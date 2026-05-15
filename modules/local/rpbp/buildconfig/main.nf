process RPBP_BUILDCONFIG {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(fasta), path(gtf)
    val(genome_name)
    val(extra_yaml)

    output:
    tuple val(meta), path("rpbp_config.yaml"), emit: config

    when:
    task.ext.when == null || task.ext.when

    script:
    def name = genome_name ?: (meta.id ?: 'reference')
    def extras = extra_yaml ?: ''
    """
    cat > rpbp_config.yaml <<EOF
genome_base_path: \$PWD/rpbp_index
genome_name: ${name}
gtf: \$PWD/${gtf}
fasta: \$PWD/${fasta}
star_index: \$PWD/rpbp_index/star
ribosomal_index: \$PWD/rpbp_index/ribosomal
ribosomal_fasta: \$PWD/${fasta}
${extras}
EOF
    """

    stub:
    """
    touch rpbp_config.yaml
    """
}
