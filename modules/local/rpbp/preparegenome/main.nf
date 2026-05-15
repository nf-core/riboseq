process RPBP_PREPAREGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(fasta), path(gtf), path(config_yaml)

    output:
    tuple val(meta), path("rpbp_index"), path(config_yaml), emit: index
    tuple val("${task.process}"), val('rpbp'), eval("prepare-rpbp-genome --version 2>&1 | sed 's/^.*rpbp //; s/ .*//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p rpbp_index
    prepare-rpbp-genome \\
        ${config_yaml} \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mkdir -p rpbp_index
    touch rpbp_index/.placeholder
    """
}
