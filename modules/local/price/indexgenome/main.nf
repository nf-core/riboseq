process PRICE_INDEXGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f293f3713b7/data' :
        'community.wave.seqera.io/library/gedi:1.0.6a--df2db6b3258d372d' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("price_index")                                                  , emit: index
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | grep -oP 'Gedi version \\K[^ ]+'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = meta.id ?: 'reference'
    """
    mkdir -p price_index
    gedi -e IndexGenome \\
        -s ${fasta} \\
        -a ${gtf} \\
        -n ${name} \\
        -f price_index \\
        -nomapping \\
        -p \\
        ${args}
    """

    stub:
    def name = meta.id ?: 'reference'
    """
    mkdir -p price_index
    touch price_index/${name}.oml
    """
}
