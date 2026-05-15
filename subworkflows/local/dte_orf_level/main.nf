//
// Tier 2 of issue #168: ORF-level DESeq2 interaction-model DTE.
//
// Per contrast, fits ~ condition + seq_type + condition:seq_type at
// ORF resolution. The Ribo-seq numerator is the per-ORF P-site matrix
// from issue #166. The RNA-seq denominator is the gene-level Salmon
// matrix, joined to each ORF via `orf_to_gene.tsv` (issue #167).
//
// Novel intergenic ORFs (no host gene) are dropped inside the R
// template because they cannot receive an RNA-seq denominator. Their
// raw counts remain available in `orf_psite_counts.tsv` for users who
// want count-only inspection.
//
// Gating happens at the call site in workflows/riboseq/main.nf:
//   --extended_orf_analysis = true  AND
//   ORF count matrix exists         AND
//   --contrasts provided            AND
//   at least one ORF caller enabled.

include { DESEQ2_ORF_DTE } from '../../../modules/local/deseq2/orf_dte'

workflow DTE_ORF_LEVEL {

    take:
    ch_contrasts          // channel: [ contrast_meta, contrast_variable, reference, target ]
    ch_orf_count_matrix   // channel: [ meta, orf_psite_counts.tsv ]    - per-ORF P-site counts (#166)
    ch_rna_counts         // channel: [ meta, gene_counts.tsv ]         - Salmon gene-level RNA-seq counts
    ch_orf_to_gene        // channel: [ meta, orf_to_gene.tsv ]         - ORF->gene mapping (#167)
    ch_samplesheet        // channel: path(sample_sheet.csv)

    main:

    ch_versions = Channel.empty()

    // Pack the per-cohort inputs into a single value channel keyed by
    // [meta, samplesheet, ribo, rna, orf_to_gene]; combined with each
    // contrast row at the call to DESEQ2_ORF_DTE below.
    ch_inputs = ch_samplesheet
        .combine(ch_orf_count_matrix.map { _meta, counts -> counts })
        .combine(ch_rna_counts.map { _meta, counts -> counts })
        .combine(ch_orf_to_gene.map { _meta, tsv -> tsv })
        .map { samplesheet, ribo, rna, o2g -> [ [id: 'allsamples'], samplesheet, ribo, rna, o2g ] }
        .first()

    DESEQ2_ORF_DTE(
        ch_contrasts,
        ch_inputs
    )
    ch_versions = ch_versions.mix(DESEQ2_ORF_DTE.out.versions)

    emit:
    results          = DESEQ2_ORF_DTE.out.results
    dispersion_plot  = DESEQ2_ORF_DTE.out.dispersion_plot
    session_info     = DESEQ2_ORF_DTE.out.session_info
    versions         = ch_versions
}
