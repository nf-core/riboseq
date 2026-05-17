//
// Run Rp-Bp's post-alignment ORF prediction chain on STAR-aligned BAMs.
//
// Chain: BUILDCONFIG -> PREPAREGENOME (once per run) -> RPBP_PREDICT_ORFS
// (per Ribo-seq sample, six rpbp tools split into one Nextflow process each).
// We deliberately avoid rpbp's `run-rpbp-pipeline` / `create-orf-profiles`
// wrappers - both invoke flexbar+bowtie+STAR for raw FASTQs, which duplicates
// the standard nf-core STAR alignment we run upstream. Splitting also gives
// independent caching per step on resume.
//
include { RPBP_BUILDCONFIG   } from '../../../modules/local/rpbp/buildconfig'
include { RPBP_PREPAREGENOME } from '../../../modules/local/rpbp/preparegenome'
include { RPBP_PREDICT_ORFS  } from '../rpbp_predict_orfs'

workflow RUN_RPBP {

    take:
    ch_bams_for_analysis  // channel: [ val(meta), path(bam), path(bai) ] - Ribo-seq only
    ch_fasta_gtf          // channel: [ val(meta), path(fasta), path(gtf) ] - first()
    extra_yaml            // val: optional extra YAML appended to the config

    main:

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

    // 3. Take the dedicated BED emissions from PREPAREGENOME (named outputs are
    //    safer than `file("${idx}/sub.bed")` under Fusion - the latter loses the
    //    s3:// scheme and the head node can't stage it).
    ch_transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed  .map { _meta, bed -> bed }.first()
    ch_orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed.map { _meta, bed -> bed }.first()
    ch_orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed  .map { _meta, bed -> bed }.first()
    ch_genome_fasta     = ch_fasta_gtf.map { _meta, fasta, _gtf -> fasta }.first()

    // 4. Per-sample ORF prediction.
    RPBP_PREDICT_ORFS (
        ch_bams_for_analysis,
        ch_transcript_bed,
        ch_orfs_genomic_bed,
        ch_orfs_exons_bed,
        ch_genome_fasta
    )

    emit:
    bed = RPBP_PREDICT_ORFS.out.predicted
}
