//
// Novel transcript discovery + hybrid GTF construction for the riboseq pipeline.
//
// Two upstream sources feed the same downstream chain
// (gffcompare -> class-code filter -> optional rRNA blacklist intersect
//  -> concat with canonical):
//   1. `--novel_gtf <file>`  : user-supplied annotation, skip StringTie.
//   2. `--skip_stringtie false`: run StringTie on RNA-seq BAMs.
//
// When neither path produces novel transcripts the caller should not invoke
// this subworkflow at all (gating lives at the call site).
//

include { BAM_STRINGTIE_MERGE        } from '../../nf-core/bam_stringtie_merge'
include { GTF_HYBRIDMERGE_GFFCOMPARE } from '../../nf-core/gtf_hybridmerge_gffcompare/main'

workflow NOVEL_TRANSCRIPT_DISCOVERY {

    take:
    ch_rnaseq_bam          // channel: [ val(meta), path(bam) ]    - genome BAMs flagged as RNA-seq
    ch_gtf                 // value channel: path(genome.gtf)      - full multi-isoform reference for gffcompare
    ch_canonical_gtf       // value channel: path(canonical.gtf)   - one-transcript-per-gene backbone for concat
    user_novel_gtf         // val: path or null                    - --novel_gtf override (skips StringTie if set)
    run_stringtie          // val: bool                            - run StringTie when no user override
    gffcompare_class_codes // val: string                          - --gffcompare_class_codes
    rrna_blacklist         // val: path or null                    - optional --rrna_blacklist for intersect
    ch_strandedness        // channel: string                      - library strandedness for strand-aware blacklist filter

    main:

    ch_novel_pre_filter = Channel.empty()

    if (user_novel_gtf) {
        ch_novel_pre_filter = Channel
            .fromPath(user_novel_gtf, checkIfExists: true)
            .combine(ch_strandedness)
            .map { gtf, strand -> [ [id: 'novel_gtf', strandedness: strand], gtf ] }
    }
    else if (run_stringtie) {
        BAM_STRINGTIE_MERGE(
            ch_rnaseq_bam,
            ch_gtf.map { gtf -> [ [:], gtf ] }
        )

        ch_novel_pre_filter = BAM_STRINGTIE_MERGE.out.stringtie_gtf
            .combine(ch_strandedness)
            .map { _meta, gtf, strand -> [ [id: 'stringtie_merge', strandedness: strand], gtf ] }
    }

    //
    // Classify novel transcripts against the full reference GTF, drop class
    // codes not in --gffcompare_class_codes, optionally subtract a strand-aware
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
        ch_gtf.map { gtf -> [ [:], gtf ] },
        ch_canonical_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] },
        Channel.value(gffcompare_class_codes),
        ch_blacklist_bed
    )

    ch_hybrid_gtf = GTF_HYBRIDMERGE_GFFCOMPARE.out.hybrid_gtf.map { _meta, gtf -> gtf }

    emit:
    hybrid_gtf = ch_hybrid_gtf // channel: path(hybrid.gtf) - canonical + filtered novel
}
