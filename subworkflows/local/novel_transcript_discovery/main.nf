//
// Novel transcript discovery + hybrid GTF construction for the riboseq pipeline.
//
// Two upstream sources feed the same downstream chain
// (gffcompare -> class-code filter -> optional rRNA blacklist intersect
//  -> concat with canonical):
//   1. `--novel_gtf <file>`  : user-supplied annotation, skip StringTie.
//   2. `--skip_stringtie false`: run StringTie on RNA-seq BAMs (preferred)
//      or fall back to Ribo-seq BAMs with tightened parameters.
//
// When neither path produces novel transcripts the caller should not invoke
// this subworkflow at all (gating lives at the call site).
//

include { BAM_STRINGTIE_MERGE        } from '../../nf-core/bam_stringtie_merge'
include { GTF_HYBRIDMERGE_GFFCOMPARE } from '../../nf-core/gtf_hybridmerge_gffcompare/main'

workflow NOVEL_TRANSCRIPT_DISCOVERY {

    take:
    ch_rnaseq_bam       // channel: [ val(meta), path(bam) ]    - genome BAMs flagged as RNA-seq
    ch_riboseq_bam      // channel: [ val(meta), path(bam) ]    - genome BAMs flagged as Ribo-seq (fallback only)
    ch_gtf              // channel: path(genome.gtf)            - full multi-isoform reference for gffcompare
    ch_canonical_gtf    // channel: path(canonical.gtf)         - one-transcript-per-gene backbone for concat
    user_novel_gtf      // val: path or null                    - --novel_gtf override (skips StringTie if set)
    run_stringtie       // val: bool                            - run StringTie when no user override
    stringtie_class_codes // val: string                        - --stringtie_class_codes
    rrna_blacklist      // val: path or null                    - optional --rrna_blacklist for intersect
    ribo_fallback_args  // val: string                          - StringTie args for Ribo-seq fallback (for warning)

    main:

    ch_novel_pre_filter = Channel.empty()

    if (user_novel_gtf) {
        ch_novel_pre_filter = Channel
            .fromPath(user_novel_gtf, checkIfExists: true)
            .map { gtf -> [ [id: 'stringtie_merge'], gtf ] }
    }
    else if (run_stringtie) {
        // Prefer RNA-seq BAMs; fall back to Ribo-seq with tightened args.
        // The fallback is signalled via meta.stringtie_fallback so the
        // `withName: STRINGTIE_STRINGTIE` block in conf/modules.config can
        // swap ext.args based on the BAM source.
        ch_rnaseq_count = ch_rnaseq_bam.count()

        ch_stringtie_rnaseq = ch_rnaseq_bam
            .map { meta, bam -> [ meta + [stringtie_fallback: false], bam ] }

        // Ribo-seq is only consumed when no RNA-seq BAMs are present.
        // `combine` with the count value channel gates the emission.
        ch_stringtie_ribo = ch_riboseq_bam
            .combine(ch_rnaseq_count)
            .filter { _meta, _bam, rna_count -> rna_count == 0 }
            .map { meta, bam, _rna_count ->
                [ meta + [stringtie_fallback: true], bam ]
            }

        ch_rnaseq_count
            .filter { it == 0 }
            .subscribe {
                log.warn "No RNA-seq BAMs available for StringTie assembly. Using Ribo-seq BAMs with the configured StringTie arguments (${ribo_fallback_args}); the pipeline default --stringtie_ribo_fallback_args is tightened to suppress P-site-pileup artefacts, override via --extra_stringtie_args at your own risk. Novel transcript assemblies may still contain artefacts; review carefully and consider supplying --novel_gtf from a dedicated RNA-seq run."
            }

        BAM_STRINGTIE_MERGE(
            ch_stringtie_rnaseq.mix(ch_stringtie_ribo),
            ch_gtf.map { gtf -> [ [:], gtf ] }
        )

        ch_novel_pre_filter = BAM_STRINGTIE_MERGE.out.stringtie_gtf
            .map { _meta, gtf -> [ [id: 'stringtie_merge'], gtf ] }
    }

    //
    // Classify novel transcripts against the full reference GTF, drop class
    // codes not in --stringtie_class_codes, optionally subtract a strand-aware
    // blacklist, then concatenate with the canonical backbone (synthesising
    // parent gene rows for any novel transcripts whose gene_id isn't already
    // in the backbone). One subworkflow call covers gffcompare + gawk filter +
    // (optional) bedtools intersect + gawk concat.
    //
    ch_blacklist_bed = rrna_blacklist
        ? Channel.fromPath(rrna_blacklist, checkIfExists: true).map { bed -> [ [id: 'rrna_blacklist'], bed ] }
        : Channel.empty()

    if (!rrna_blacklist) {
        log.info "No rRNA/repeat blacklist supplied via --rrna_blacklist; skipping post-assembly blacklist intersect."
    }

    GTF_HYBRIDMERGE_GFFCOMPARE(
        ch_novel_pre_filter,
        ch_gtf.map { gtf -> [ [:], gtf ] }.first(),
        ch_canonical_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] }.first(),
        Channel.value(stringtie_class_codes),
        ch_blacklist_bed
    )

    ch_hybrid_gtf = GTF_HYBRIDMERGE_GFFCOMPARE.out.hybrid_gtf.map { _meta, gtf -> gtf }

    emit:
    hybrid_gtf = ch_hybrid_gtf // channel: path(hybrid.gtf) - canonical + filtered novel
}
