process ORF_CATALOGUE_EXTRACTAA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e4a7dcecd7abc2ee2e0af0c3b11d507f9f2ac523ebe84cdc42fda3f3157c986/data' :
        'community.wave.seqera.io/library/orf_catalogue_merge:8934e10959c9dfa3' }"

    input:
    tuple val(meta), path(catalogue_bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("orf_catalogue.faa"), emit: faa
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Index FASTA if no .fai alongside (pysam will fail otherwise).
    if [ ! -e ${fasta}.fai ]; then
        python -c "import pysam; pysam.faidx('${fasta}')"
    fi
    orf_catalogue_extractaa.py \\
        --bed12 ${catalogue_bed} \\
        --fasta ${fasta} \\
        --out-faa orf_catalogue.faa
    """

    stub:
    """
    touch orf_catalogue.faa
    """
}
