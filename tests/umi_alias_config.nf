include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_HYBRID     } from '../subworkflows/nf-core/bam_dedup_umi'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_HYBRID_RNA } from '../subworkflows/nf-core/bam_dedup_umi'

workflow UMI_ALIAS_CONFIG {

    take:
    ch_bam
    ch_fasta_fai
    umi_dedup_tool
    umitools_dedup_stats
    ch_transcriptome_bam
    ch_transcript_fasta_fai
    umitools_dedup_primary_only

    main:
    BAM_DEDUP_UMI_HYBRID(
        ch_bam,
        ch_fasta_fai,
        umi_dedup_tool,
        umitools_dedup_stats,
        ch_transcriptome_bam,
        ch_transcript_fasta_fai,
        umitools_dedup_primary_only
    )

    BAM_DEDUP_UMI_HYBRID_RNA(
        ch_bam,
        ch_fasta_fai,
        umi_dedup_tool,
        umitools_dedup_stats,
        ch_transcriptome_bam,
        ch_transcript_fasta_fai,
        umitools_dedup_primary_only
    )

    emit:
    hybrid_transcriptome_bam          = BAM_DEDUP_UMI_HYBRID.out.transcriptome_bam
    hybrid_transcriptome_dedup_bam    = BAM_DEDUP_UMI_HYBRID.out.transcriptome_dedup_bam
    hybrid_transcriptome_sorted_bam   = BAM_DEDUP_UMI_HYBRID.out.transcriptome_sorted_bam
    hybrid_transcriptome_filtered_bam = BAM_DEDUP_UMI_HYBRID.out.transcriptome_filtered_bam
    rna_transcriptome_bam             = BAM_DEDUP_UMI_HYBRID_RNA.out.transcriptome_bam
    rna_transcriptome_dedup_bam       = BAM_DEDUP_UMI_HYBRID_RNA.out.transcriptome_dedup_bam
    rna_transcriptome_sorted_bam      = BAM_DEDUP_UMI_HYBRID_RNA.out.transcriptome_sorted_bam
    rna_transcriptome_filtered_bam    = BAM_DEDUP_UMI_HYBRID_RNA.out.transcriptome_filtered_bam
}
