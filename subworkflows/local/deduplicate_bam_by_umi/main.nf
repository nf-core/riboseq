include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'

workflow DEDUPLICATE_BAM_BY_UMI {

    // Deduplicate genome BAM file before downstream analysis
    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
        ch_genome_bam.join(ch_genome_bam_index, by: [0]),
        params.umitools_dedup_stats
    )
    ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
    ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
    ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
    ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
    ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
    if (params.bam_csi_index) {
        ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
    }
    ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)

    // Co-ordinate sort, index and run stats on transcriptome BAM
    BAM_SORT_STATS_SAMTOOLS (
        ch_transcriptome_bam,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] }
    )
    ch_transcriptome_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_transcriptome_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

    // Deduplicate transcriptome BAM file before read counting with Salmon
    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
        ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
        params.umitools_dedup_stats
    )

    // Name sort BAM before passing to Salmon
    SAMTOOLS_SORT (
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
    )

    // Only run prepare_for_rsem.py on paired-end BAM files
    SAMTOOLS_SORT
        .out
        .bam
        .branch {
            meta, bam ->
        single_end: meta.single_end
            return [ meta, bam ]
        paired_end: !meta.single_end
            return [ meta, bam ]
        }
        .set { ch_umitools_dedup_bam }

    // Fix paired-end reads in name sorted BAM file
    // See: https://github.com/nf-core/rnaseq/issues/828
    UMITOOLS_PREPAREFORSALMON (
        ch_umitools_dedup_bam.paired_end
    )
    ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

    ch_umitools_dedup_bam
        .single_end
        .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
        .set { ch_transcriptome_bam }
}
