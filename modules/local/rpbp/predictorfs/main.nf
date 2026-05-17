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
    # rpbp 4.x with `keep_riboseq_multimappers: True` looks for
    # without-rrna-mapping/<sample>.bam (filenames.py get_riboseq_bam).
    mkdir -p rpbp_out without-rrna-mapping
    ln -sf ../${bam} without-rrna-mapping/${prefix}.bam
    ln -sf ../${bai} without-rrna-mapping/${prefix}.bam.bai

    cp ${config_yaml} predict_config.yaml
    cat >> predict_config.yaml <<EOF
riboseq_data: .
riboseq_samples:
  ${prefix}: ${prefix}.bam
keep_riboseq_multimappers: True
EOF

    run-rpbp-pipeline \\
        ${bam} \\
        predict_config.yaml \\
        ${prefix} \\
        --num-cpus ${task.cpus} \\
        ${args}

    find . \\( -name "*predicted-orfs*.bed.gz" -o -name "*predicted-orfs*.tab.gz" \\) \\
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
