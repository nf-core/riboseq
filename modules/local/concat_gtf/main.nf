process CONCAT_GTF {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(canonical_gtf), path(novel_gtf, stageAs: 'novel.gtf')

    output:
    tuple val(meta), path("*.gtf")    , emit: gtf
    tuple val("${task.process}"), val('coreutils'), eval("sort --version | head -n1 | sed 's/.*) //; s/ .*//'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def novel  = novel_gtf ? "novel.gtf" : ""
    """
    cat ${canonical_gtf} ${novel} \\
        | grep -v '^#' \\
        | sort -k1,1 -k4,4n \\
        > ${prefix}.gtf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    """
}
