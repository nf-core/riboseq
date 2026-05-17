process RPBP_ESTIMATEMETAGENEBAYESFACTORS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(metagene_profile)

    output:
    tuple val(meta), path("${prefix}.metagene-periodicity-bayes-factors.csv.gz"), emit: bayes_factors
    tuple val("${task.process}"), val('rpbp'), eval("estimate-metagene-profile-bayes-factors --version 2>&1 | sed 's/^.*rpbp //; s/ .*//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    estimate-metagene-profile-bayes-factors \\
        ${metagene_profile} \\
        ${prefix}.metagene-periodicity-bayes-factors.csv.gz \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.metagene-periodicity-bayes-factors.csv.gz
    """
}
