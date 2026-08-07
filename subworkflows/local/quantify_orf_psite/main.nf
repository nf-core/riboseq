//
// ORF-level in-frame P-site quantification.
//
// Chains:
//   1. CUSTOM_BED12CODONPOSITIONS: expand the BED12 catalogue into codon-start BED6.
//   2. QUANTIFY_INFRAME_PSITE_PLASTID: per-sample bedtools intersect + groupby
//      against plastid wiggle tracks (shared with the gene-level path).
//   3. ORF_COUNT_MATRIX: assemble per-sample TSVs into one ORF x sample
//      count matrix, zero-filled for ORFs absent from a sample.
//
// Gating is done at the call site in workflows/riboseq/main.nf:
//   --extended_orf_analysis = true  AND  a caller ran  AND  plastid is enabled.
//
// The matrix is the input for the downstream DTE step. This
// subworkflow does not run DTE itself.

include { CUSTOM_BED12CODONPOSITIONS     } from '../../../modules/nf-core/custom/bed12codonpositions/main'
include { QUANTIFY_INFRAME_PSITE_PLASTID } from '../../../modules/local/quantify_inframe_psite_plastid'
include { ORF_COUNT_MATRIX               } from '../../../modules/local/orf_count_matrix'

workflow QUANTIFY_ORF_PSITE {

    take:
    ch_catalogue_bed   // channel: [ val(meta), path(orf_catalogue.bed12) ]
    ch_psite_tracks    // channel: [ val(meta), path(forward.bedgraph), path(reverse.bedgraph) ] - riboseq samples only

    main:

    ch_versions = Channel.empty()

    // 1. Expand catalogue to codon-start BED6 (cohort-level, runs once).
    CUSTOM_BED12CODONPOSITIONS ( ch_catalogue_bed )
    ch_versions = ch_versions.mix(CUSTOM_BED12CODONPOSITIONS.out.versions)

    // `feature: 'orf'` sets the output-filename suffix in the shared counter.
    // Must be a value channel: it is paired with the per-sample track channel
    // below, and a single-item queue would serve only one sample.
    ch_inframe_psites = CUSTOM_BED12CODONPOSITIONS.out.bed
        .map { meta, bed -> [ meta + [feature: 'orf'], bed ] }
        .first()

    // 2. Per-sample P-site counting (shared module with the gene-level path).
    QUANTIFY_INFRAME_PSITE_PLASTID ( ch_psite_tracks, ch_inframe_psites )
    ch_versions = ch_versions.mix(QUANTIFY_INFRAME_PSITE_PLASTID.out.versions)

    // 3. Collect all per-sample TSVs into one list, then pair with the
    //    cohort-level catalogue BED on a synthetic key so the pairing does not
    //    depend on emission order.
    ch_tsvs_collected = QUANTIFY_INFRAME_PSITE_PLASTID.out.counts
        .map { _meta, tsv -> tsv }
        .collect()
        .map { tsvs -> [ 'allsamples', tsvs ] }

    ch_bed_keyed = ch_catalogue_bed
        .map { _meta, bed -> [ 'allsamples', bed ] }

    ch_matrix_in = ch_tsvs_collected
        .combine( ch_bed_keyed, by: 0 )
        .map { _key, tsvs, bed -> [ [ id: 'allsamples' ], tsvs, bed ] }

    ORF_COUNT_MATRIX ( ch_matrix_in )
    ch_versions = ch_versions.mix(ORF_COUNT_MATRIX.out.versions)

    emit:
    orf_count_matrix = ORF_COUNT_MATRIX.out.matrix             // [ meta, orf_psite_counts.tsv ]
    inframe_psites   = CUSTOM_BED12CODONPOSITIONS.out.bed      // [ meta, *.bed ]
    versions         = ch_versions
}
