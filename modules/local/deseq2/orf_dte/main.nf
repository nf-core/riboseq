// Tier 2 of issue #168: per-ORF DESeq2 interaction-model DTE.
//
// Fits ~ condition + seq_type + condition:seq_type per ORF using the
// per-ORF P-site count matrix as the Ribo-seq numerator and the
// gene-level RNA-seq counts as the denominator (joined to each ORF via
// orf_to_gene.tsv; novel intergenic ORFs with no host gene are dropped).
//
// The container is reused from the existing deltaTE module: it already
// carries DESeq2 + apeglm + data.table + dplyr + readr + tibble. No new
// container build is required.
//
// Row-independence caveat: when multiple ORFs share the same gene-level
// RNA-seq denominator row, those rows are not statistically independent
// after the gene-level RNA counts are joined in. We proceed per-ORF
// with this assumption explicitly documented (see docs/usage.md). A
// statistically rigorous solution (Fishpond/Swish with inferential
// replicates) is out of scope.

process DESEQ2_ORF_DTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b39a67d1085303bdc1a56f1ab0da64673269e065297d2c58f66a42a025c97e84/data':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-apeglm_bioconductor-complexheatmap_r-data.table_pruned:fdc920cced486331' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(ribo_counts), path(rna_counts), path(orf_to_gene)

    output:
    tuple val(meta), path("*.orf_dte.results.tsv")          , emit: results
    tuple val(meta), path("*.orf_dte.dispersions.png")      , emit: dispersion_plot, optional: true
    tuple val(meta), path("*.orf_dte.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_orf_dte.R'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.orf_dte.results.tsv
    touch ${prefix}.orf_dte.dispersions.png
    touch ${prefix}.orf_dte.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: 1.42.0
    END_VERSIONS
    """
}
