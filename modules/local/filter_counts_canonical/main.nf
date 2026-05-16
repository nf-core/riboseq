process FILTER_COUNTS_CANONICAL {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(counts)
    tuple val(meta2), path(canonical_gtf)

    output:
    tuple val(meta), path("*.canonical_filtered.tsv"), emit: counts
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id ?: counts.baseName}"
    """
    awk -F'\\t' '
        BEGIN { OFS=FS }
        FNR == NR {
            if (\$0 ~ /^#/) { next }
            if (\$3 != "gene" && \$3 != "transcript") { next }
            if (match(\$9, /gene_id "([^"]+)"/, m)) keep[m[1]] = 1
            next
        }
        FNR == 1 { print; next }
        (\$1 in keep) { print }
    ' "${canonical_gtf}" "${counts}" > "${prefix}.canonical_filtered.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion 2>/dev/null | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id ?: counts.baseName}"
    """
    touch "${prefix}.canonical_filtered.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion 2>/dev/null | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
