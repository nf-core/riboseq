//
// Equalise RNA-seq read lengths to match Ribo-seq read lengths
//

include { SEQKIT_STATS                      } from '../../../modules/nf-core/seqkit/stats'
include { TRIMGALORE as TRIMGALORE_HARDTRIM } from '../../../modules/nf-core/trimgalore'

// Helper function to normalize meta fields that may be null, empty list, or value
// Also converts string numbers to integers for trim_length compatibility
def normalizeMetaField(value) {
    def result = value instanceof List ? (value ? value[0] : null) : value
    // Convert string numbers to integers (for trim_length from samplesheet)
    if (result instanceof String && result.isInteger()) {
        return result.toInteger()
    }
    return result
}

workflow FASTQ_EQUALISE_READ_LENGTHS {
    take:
    ch_reads                        // channel: [ val(meta), [ path(reads) ] ]
    equalise_read_lengths_target    // integer: Global target length (optional fallback)

    main:
    ch_versions = Channel.empty()

    // Branch reads by sample type
    ch_reads
        .branch { meta, reads ->
            riboseq: meta.sample_type == 'riboseq'
                return [ meta, reads ]
            rnaseq: meta.sample_type == 'rnaseq'
                return [ meta, reads ]
            other: true
                return [ meta, reads ]
        }
        .set { ch_reads_by_type }

    // Determine trim length for each rnaseq sample based on priority:
    // 1. Per-sample samplesheet trim_length on the RNA-seq sample itself
    // 2. Global equalise_read_lengths_target parameter
    // 3. Derived from paired Ribo-seq via seqkit stats (only if needed)

    // First pass: resolve samples that have explicit trim_length or can use global target
    // Note: trim_length may be null, empty list [], or an integer depending on samplesheet parsing
    ch_rnaseq_first_pass = ch_reads_by_type.rnaseq
        .map { meta, reads ->
            def trim_len = normalizeMetaField(meta.trim_length)
            if (trim_len) {
                // Priority 1: Per-sample trim_length on RNA-seq sample
                return [ meta + [resolved_trim_length: trim_len], reads ]
            } else if (equalise_read_lengths_target) {
                // Priority 2: Global target length
                return [ meta + [resolved_trim_length: equalise_read_lengths_target], reads ]
            }
            // Needs to be derived from paired riboseq
            return [ meta, reads ]
        }
        .branch { meta, reads ->
            has_length: meta.resolved_trim_length != null
                return [ meta, reads ]
            needs_pairing: true
                return [ meta, reads ]
        }

    // Only run seqkit stats if there are samples that need pairing
    // Filter riboseq samples to only those with pair values that match rnaseq samples needing pairing
    ch_rnaseq_pairs_needed = ch_rnaseq_first_pass.needs_pairing
        .filter { meta, reads -> meta.pair != null }
        .map { meta, reads -> meta.pair }
        .collect()
        .map { pairs_list -> [pairs_list] }  // Wrap in list to prevent spreading in combine
        .ifEmpty([[]])

    ch_riboseq_to_stats = ch_reads_by_type.riboseq
        .combine(ch_rnaseq_pairs_needed)
        .filter { meta, reads, pairs_needed ->
            meta.pair != null && pairs_needed.contains(meta.pair)
        }
        .map { meta, reads, pairs_needed -> [ meta, reads ] }

    // Run seqkit stats only on riboseq samples that have matching rnaseq pairs needing length
    SEQKIT_STATS(ch_riboseq_to_stats)

    // Parse seqkit stats output to extract average length
    ch_riboseq_lengths = SEQKIT_STATS.out.stats
        .map { meta, stats_file ->
            def lines = stats_file.text.readLines()
            def meta_trim_len = normalizeMetaField(meta.trim_length)
            if (lines.size() > 1) {
                def header = lines[0].split('\t')
                def values = lines[1].split('\t')
                def avg_len_idx = header.findIndexOf { it == 'avg_len' }
                def avg_len = avg_len_idx >= 0 ? Math.round(values[avg_len_idx].toFloat()) : null
                return [ meta.pair, meta_trim_len ?: avg_len ]
            }
            return [ meta.pair, null ]
        }

    // Join rnaseq samples needing pairing with their riboseq-derived lengths
    ch_rnaseq_from_pairing = ch_rnaseq_first_pass.needs_pairing
        .map { meta, reads -> [ meta.pair ?: 'NO_PAIR', meta, reads ] }
        .join(ch_riboseq_lengths, by: 0, remainder: true)
        .map { pair, meta, reads, riboseq_trim_len ->
            if (riboseq_trim_len == null) {
                error "RNA-seq sample '${meta.id}' (pair: ${meta.pair ?: 'none'}) has no matching Ribo-seq sample and no --equalise_read_lengths_target specified"
            }
            return [ meta + [resolved_trim_length: riboseq_trim_len], reads ]
        }

    // Combine all rnaseq samples with resolved trim lengths
    ch_rnaseq_to_trim = ch_rnaseq_first_pass.has_length
        .mix(ch_rnaseq_from_pairing)

    // For paired-end, take only R1 and mark as single-end
    ch_rnaseq_for_hardtrim = ch_rnaseq_to_trim
        .map { meta, reads ->
            def new_meta = meta + [
                single_end: true
            ]
            // For paired-end, take only R1
            def new_reads = meta.single_end ? reads : reads[0]
            return [ new_meta, new_reads ]
        }

    // Run trimgalore hardtrim on rnaseq samples
    TRIMGALORE_HARDTRIM(ch_rnaseq_for_hardtrim)

    // Collect versions - use ifEmpty to handle processes that may not run
    ch_versions = ch_versions
        .mix(SEQKIT_STATS.out.versions.first().ifEmpty([]))
        .mix(TRIMGALORE_HARDTRIM.out.versions.first().ifEmpty([]))
        .filter { it }  // Remove empty elements

    // Combine trimmed rnaseq with riboseq and other samples (normalizing meta)
    // Normalize trim_length in meta (convert empty list to null), but only if it exists
    ch_reads_out = TRIMGALORE_HARDTRIM.out.reads
        .mix(ch_reads_by_type.riboseq.map { meta, reads ->
            if (meta.containsKey('trim_length')) {
                return [ meta + [trim_length: normalizeMetaField(meta.trim_length)], reads ]
            }
            return [ meta, reads ]
        })
        .mix(ch_reads_by_type.other.map { meta, reads ->
            if (meta.containsKey('trim_length')) {
                return [ meta + [trim_length: normalizeMetaField(meta.trim_length)], reads ]
            }
            return [ meta, reads ]
        })

    emit:
    reads           = ch_reads_out                          // channel: [ val(meta), [ path(reads) ] ]
    riboseq_stats   = SEQKIT_STATS.out.stats                // channel: [ val(meta), path(stats) ] - empty if all samples have explicit trim_length
    trim_log        = TRIMGALORE_HARDTRIM.out.log           // channel: [ val(meta), path(log) ] - may be empty in hardtrim mode
    versions        = ch_versions                           // channel: [ versions.yml ]
}
