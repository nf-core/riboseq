process ORF_NORMALISE_RIBOTISH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5dca02d2332863bdfef5216c395c373fb8c11e1c3de090b0b9f2f9b7a19afaa6/data' :
        'community.wave.seqera.io/library/orf_catalogue_normalise:ce5d6a5f27e4a6e4' }"

    input:
    tuple val(meta), path(pred_txt)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.normalised.bed12"), emit: bed12
    tuple val(meta), path("*.normalised.tsv")  , emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = meta.id ?: 'unknown'
    """
    orf_normalise_ribotish.py \\
        --input ${pred_txt} \\
        --gtf ${gtf} \\
        --sample-id ${sample_id} \\
        --out-bed ${prefix}.normalised.bed12 \\
        --out-tsv ${prefix}.normalised.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.normalised.bed12
    touch ${prefix}.normalised.tsv
    """
}
