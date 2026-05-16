process CONCAT_GTF {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(canonical_gtf), path(novel_gtf, stageAs: 'novel.gtf')

    output:
    tuple val(meta), path("*.gtf")    , emit: gtf
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion 2>/dev/null | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def novel  = novel_gtf ? "novel.gtf" : ""
    """
    cat ${canonical_gtf} ${novel} | grep -v '^#' > combined.gtf

    awk -F'\\t' 'BEGIN { OFS=FS }
        \$3 == "gene" {
            if (match(\$9, /gene_id "([^"]+)"/, m)) gene_seen[m[1]] = 1
            print; next
        }
        \$3 == "transcript" {
            if (match(\$9, /gene_id "([^"]+)"/, m)) {
                gid = m[1]
                if (!(gid in gene_seen) && !(gid in synth)) {
                    print \$1, \$2, "gene", \$4, \$5, ".", \$7, ".", "gene_id \\"" gid "\\";"
                    synth[gid] = 1
                }
            }
            print; next
        }
        { print }
    ' combined.gtf | sort -k1,1 -k4,4n > ${prefix}.gtf

    rm combined.gtf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    """
}
