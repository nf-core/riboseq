process ORF_NORMALISE {
    tag "${meta.id}.${meta.caller}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a17ff642fb8c74fbf9feece8df823d78bcecca69f76770aa585e7468e6a9187/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:3ce680b21acf323f' }"

    input:
    tuple val(meta), path(orfs_table)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.normalised.bed12"), emit: bed12
    tuple val(meta), path("*.normalised.tsv")  , emit: tsv
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    sample_id = meta.id ?: 'unknown'
    caller = meta.caller
    template 'orf_normalise.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.normalised.bed12
    touch ${prefix}.normalised.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
