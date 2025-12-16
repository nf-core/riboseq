process DESEQ2_DELTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b39a67d1085303bdc1a56f1ab0da64673269e065297d2c58f66a42a025c97e84/data':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-apeglm_bioconductor-complexheatmap_r-data.table_pruned:fdc920cced486331' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.translation.deltate.results.tsv")   , emit: translation
    tuple val(meta), path("*.translated_mRNA.deltate.results.tsv"), emit: translated_mrna
    tuple val(meta), path("*.total_mRNA.deltate.results.tsv")    , emit: total_mrna
    tuple val(meta), path("*.dtegs.deltate.genes.tsv")           , emit: dtegs
    tuple val(meta), path("*.mRNA_abundance.deltate.genes.tsv")  , emit: mrna_abundance
    tuple val(meta), path("*.translation.deltate.genes.tsv")     , emit: translation_genes
    tuple val(meta), path("*.intensified.deltate.genes.tsv")     , emit: intensified
    tuple val(meta), path("*.buffering.deltate.genes.tsv")       , emit: buffering
    tuple val(meta), path("*.fold_change.png")                  , emit: fold_change_plot
    tuple val(meta), path("*.interaction_p_distribution.png")   , emit: interaction_p_distribution_plot, optional: true
    tuple val(meta), path("*.residual_distribution_summary.png") , emit: residual_distribution_summary_plot, optional: true
    tuple val(meta), path("*.residual_vs_fitted.png")           , emit: residual_vs_fitted_plot, optional: true
    tuple val(meta), path("*.effect_size_rna_vs_ribo.png")      , emit: effect_size_rna_vs_ribo_plot, optional: true
    tuple val(meta), path("*.effect_size_volcano.png")          , emit: effect_size_volcano_plot, optional: true
    tuple val(meta), path("*.pca_ribo.png"), path("*.pca_ribo.tsv")                                          , emit: pca_ribo , optional: true
    tuple val(meta), path("*.pca_rna.png"), path("*.pca_rna.tsv")                                            , emit: pca_rna  , optional: true
    tuple val(meta), path("*.heatmap.png"), path("*.heatmap_zscores.tsv"), path("*.heatmap_annotations.tsv") , emit: heatmap  , optional: true
    tuple val(meta), path("*.DESeqDataSet.rds")                 , emit: rdata
    tuple val(meta), path("*.R_sessionInfo.log")                , emit: session_info
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_deltate.R'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.translation.deltate.results.tsv
    touch ${prefix}.translated_mRNA.deltate.results.tsv
    touch ${prefix}.total_mRNA.deltate.results.tsv
    touch ${prefix}.dtegs.deltate.genes.tsv
    touch ${prefix}.mRNA_abundance.deltate.genes.tsv
    touch ${prefix}.translation.deltate.genes.tsv
    touch ${prefix}.intensified.deltate.genes.tsv
    touch ${prefix}.buffering.deltate.genes.tsv
    touch ${prefix}.fold_change.png
    touch ${prefix}.interaction_p_distribution.png
    touch ${prefix}.residual_distribution_summary.png
    touch ${prefix}.residual_vs_fitted.png
    touch ${prefix}.effect_size_rna_vs_ribo.png
    touch ${prefix}.effect_size_volcano.png
    touch ${prefix}.pca_ribo.png
    touch ${prefix}.pca_ribo.tsv
    touch ${prefix}.pca_rna.png
    touch ${prefix}.pca_rna.tsv
    touch ${prefix}.heatmap.png
    touch ${prefix}.heatmap_zscores.tsv
    touch ${prefix}.heatmap_annotations.tsv
    touch ${prefix}.DESeqDataSet.rds
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: 1.42.0
        r-ggplot2: 3.4.4
    END_VERSIONS
    """
}
