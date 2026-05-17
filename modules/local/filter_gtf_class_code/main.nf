process FILTER_GTF_CLASS_CODE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(annotated_gtf)
    val   class_codes

    output:
    tuple val(meta), path("*.gtf")    , emit: gtf
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion 2>/dev/null | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def codes_csv    = (class_codes instanceof List) ? class_codes.join(',') : "${class_codes}"
    """
    # Single-pass: requires transcript rows before their child exons (gffcompare order).
    awk -v codes="${codes_csv}" 'BEGIN {
        n = split(codes, arr, /,[[:space:]]*/);
        for (i = 1; i <= n; i++) keep_code[arr[i]] = 1;
    }
    /^#/ { next }
    \$7 == "." { next }
    {
        if (match(\$0, /transcript_id "([^"]+)"/, t)) {
            if (match(\$0, /class_code "([^"]+)"/, c)) {
                if (c[1] in keep_code) keep_tx[t[1]] = 1;
            }
            if (t[1] in keep_tx) print \$0;
        }
    }' ${annotated_gtf} > ${prefix}.gtf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    """
}
