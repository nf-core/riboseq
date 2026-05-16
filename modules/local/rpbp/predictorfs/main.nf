process RPBP_PREDICTORFS {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(bam), path(bai), path(prepared_index), path(config_yaml)

    output:
    tuple val(meta), path("*.predicted-orfs.bed.gz"), emit: bed, optional: true
    tuple val(meta), path("*.predicted-orfs.tab.gz"), emit: tab, optional: true
    tuple val(meta), path("rpbp_out")               , emit: outdir
    tuple val("${task.process}"), val('rpbp'), eval("predict-translated-orfs --version 2>&1 | sed 's/^.*rpbp //; s/ .*//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rpbp_out

    # predict-translated-orfs requires a `riboseq_data` mapping in the config
    # so it can resolve <sample-name> -> BAM. BUILDCONFIG renders the shared
    # genome config but doesn't know the sample list, so we append the
    # current sample's entry here per-task.
    cp ${config_yaml} predict_config.yaml
    cat >> predict_config.yaml <<EOF
riboseq_data:
  ${prefix}: ${bam}
EOF

    predict-translated-orfs \\
        predict_config.yaml \\
        ${prefix} \\
        --num-cpus ${task.cpus} \\
        ${args}

    # Surface the per-sample predicted-orfs outputs at top-level for easy publishing.
    find rpbp_out \\( -name "*predicted-orfs*.bed.gz" -o -name "*predicted-orfs*.tab.gz" \\) \\
        -exec cp -L {} . \\; || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rpbp_out
    touch ${prefix}.predicted-orfs.bed.gz
    touch ${prefix}.predicted-orfs.tab.gz
    """
}
