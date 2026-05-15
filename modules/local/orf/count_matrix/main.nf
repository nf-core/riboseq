// Merge per-sample ORF P-site count TSVs into a single ORF x sample
// matrix. The matrix has one row per ORF in the catalogue (zero-filled
// for samples with no observed P-sites) and one column per sample, in
// lexicographic sample-id order.

process ORF_COUNT_MATRIX {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(per_sample_counts, stageAs: 'counts/*'), path(orf_catalogue_bed12)

    output:
    tuple val(meta), path("orf_psite_counts.tsv"), emit: matrix
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_orf_count_matrix.py \\
        counts/*.tsv \\
        --orf-list $orf_catalogue_bed12 \\
        -o orf_psite_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    touch orf_psite_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """
}
