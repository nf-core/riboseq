process PRICE_INDEXGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/436a8b9b8bf7be4bb2122a21daef20c11b5b01353e1232ced2702f87792ff71e/data' :
        'community.wave.seqera.io/library/price_indexgenome:21e4498061541eb9' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("price_index")                                                  , emit: index
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | sed -n 's/.*Gedi version \\([^ ]*\\).*/\\1/p' | head -n 1"), emit: versions, topic: versions

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
