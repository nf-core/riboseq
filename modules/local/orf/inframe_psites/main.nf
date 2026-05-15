// Convert a BED12 ORF catalogue to a BED6 of in-frame codon-start
// (P-site) positions. Frame is defined by the ORF's own start codon,
// independent of any GTF `phase` annotation (catalogue ORFs may originate
// from novel transcripts where phase is undefined). The output is consumed
// by QUANTIFY_INFRAME_PSITE_ORF.

process ORF_INFRAME_PSITES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(catalogue_bed12)

    output:
    tuple val(meta), path("*.orf_inframe_psites.bed"), emit: bed
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    orf_bed12_to_inframe_psites.py \\
        $catalogue_bed12 \\
        -o ${prefix}.orf_inframe_psites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orf_inframe_psites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """
}
