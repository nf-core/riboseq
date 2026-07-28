process GEDI_INDEXGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/bae29fa913dea79a3dcdbfbf544f0391f82bbfdbf3e6430f71db45ba21d6cf79/data' :
        'community.wave.seqera.io/library/gedi_indexgenome:cfca16738f306c86' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: index
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | sed -n 's/.*Gedi version \\([^ ]*\\).*/\\1/p' | head -n 1"), topic: versions, emit: versions_gedi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = meta.id ?: 'reference'
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    # GEDI records absolute paths inside the .oml and .fi it writes, including one
    # pointing at the FASTA given to -s. A staged input sits at the task work
    # directory root, which is not part of the emitted index directory, so on
    # executors where work directories are not a shared filesystem (AWS Batch with
    # Fusion, and cloud executors generally) that path is never materialised and
    # downstream PRICE dies in PriceCodonInference reading sequence. Index from a
    # copy inside the index directory so the recorded path travels with the index.
    cp -L ${fasta} ${prefix}/${fasta.name}

    gedi -e IndexGenome \\
        -s ${prefix}/${fasta.name} \\
        -a ${gtf} \\
        -n ${name} \\
        -f ${prefix} \\
        -o ${prefix}/${name}.oml \\
        -nomapping \\
        -p \\
        ${args}
    """

    stub:
    def name = meta.id ?: 'reference'
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${name}.oml
    """
}
