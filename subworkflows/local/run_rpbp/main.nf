//
// Run Rp-Bp end-to-end for the Ribo-seq samples.
//
// Chain: BUILD_RPBP_CONFIG -> RPBP_PREPAREGENOME -> RPBP_PREDICTORFS
// (per sample). The same YAML config and prepared index are reused
// across samples.
//
include { RPBP_BUILDCONFIG   } from '../../../modules/local/rpbp/buildconfig'
include { RPBP_PREPAREGENOME } from '../../../modules/local/rpbp/preparegenome'
include { RPBP_PREDICTORFS   } from '../../../modules/local/rpbp/predictorfs'

workflow RUN_RPBP {

    take:
    ch_bams_for_analysis  // channel: [ val(meta), path(bam), path(bai) ] - Ribo-seq only
    ch_fasta_gtf          // channel: [ val(meta), path(fasta), path(gtf) ] - first()
    extra_yaml            // val: optional extra YAML appended to the config

    main:

    ch_versions = Channel.empty()

    // 1. Build the Rp-Bp YAML config from pipeline inputs.
    RPBP_BUILDCONFIG(
        ch_fasta_gtf,
        'reference',
        extra_yaml
    )

    // 2. Prepare the Rp-Bp index once per pipeline invocation.
    ch_preparegenome_in = ch_fasta_gtf
        .combine(RPBP_BUILDCONFIG.out.config.map { _meta, cfg -> cfg })
        .map { meta, fasta, gtf, cfg -> [ meta, fasta, gtf, cfg ] }

    RPBP_PREPAREGENOME(ch_preparegenome_in)

    // 3. Per-sample ORF prediction. Combine each sample BAM with the shared
    //    prepared index + config (single value channel, broadcast to all samples).
    ch_index = RPBP_PREPAREGENOME.out.index
        .map { _meta, idx, cfg -> [ idx, cfg ] }
        .first()

    ch_predict_in = ch_bams_for_analysis
        .combine(ch_index)
        .map { meta, bam, bai, idx, cfg -> [ meta, bam, bai, idx, cfg ] }

    RPBP_PREDICTORFS(ch_predict_in)

    emit:
    bed      = RPBP_PREDICTORFS.out.bed
    tab      = RPBP_PREDICTORFS.out.tab
    outdir   = RPBP_PREDICTORFS.out.outdir
    versions = ch_versions
}
