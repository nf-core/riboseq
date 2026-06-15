process CUSTOM_ORFCOLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15f690593e9d2bd7ebeb72fff4068eeae9cdfc909103de457b81a25d29cfbef3/data' :
        'community.wave.seqera.io/library/python_pyyaml:25b0f22c7e7bf847' }"

    input:
    tuple val(meta), path(bed12, stageAs: 'input/*'), path(catalogue_tsv, stageAs: 'input/*'), path(orf_to_gene_tsv, stageAs: 'input/*'), path(mqc_tsv, stageAs: 'input/*'), path(aa_fasta, stageAs: 'input/*'), path(cluster_tsv, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.catalogue.bed12")   , emit: bed12
    tuple val(meta), path("${prefix}.catalogue.tsv")     , emit: catalogue_tsv
    tuple val(meta), path("${prefix}.orf_to_gene.tsv")   , emit: orf_to_gene_tsv
    tuple val(meta), path("${prefix}.catalogue.mqc.tsv") , emit: multiqc
    tuple val(meta), path("${prefix}.catalogue.aa.fasta"), emit: aa_fasta
    path "versions.yml"                                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''
    template 'orfcollapse.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.catalogue.bed12 ${prefix}.catalogue.tsv ${prefix}.orf_to_gene.tsv ${prefix}.catalogue.mqc.tsv ${prefix}.catalogue.aa.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
