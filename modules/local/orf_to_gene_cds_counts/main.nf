// Tier 1 of issue #168: re-aggregate per-ORF P-site counts to gene
// level using ONLY canonical_cds ORFs. The output is long-format
// `sample<TAB>gene_id<TAB>count`, matching the shape of the per-sample
// plastid `gene_inframe_psite_counts.tsv` files, so the same
// REPLACE_RIBOSEQ_COUNTS_IN_MATRIX awk substitution step can consume
// it without changes. uORF / dORF / novel_u / smORF / other ORFs are
// excluded from the sum to keep the translational numerator clean.

process ORF_TO_GENE_CDS_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(orf_count_matrix), path(orf_to_gene_tsv), path(orf_catalogue_tsv)

    output:
    tuple val(meta), path("gene_cds_psite_counts.tsv"), emit: gene_counts
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    orf_to_gene_cds_counts.py \\
        --orf-counts ${orf_count_matrix} \\
        --orf-to-gene ${orf_to_gene_tsv} \\
        --orf-catalogue ${orf_catalogue_tsv} \\
        -o gene_cds_psite_counts.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    touch gene_cds_psite_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """
}
