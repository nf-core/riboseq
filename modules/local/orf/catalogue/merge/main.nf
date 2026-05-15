process ORF_CATALOGUE_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e4a7dcecd7abc2ee2e0af0c3b11d507f9f2ac523ebe84cdc42fda3f3157c986/data' :
        'community.wave.seqera.io/library/orf_catalogue_merge:8934e10959c9dfa3' }"

    input:
    tuple val(meta), path(bed12s, stageAs: 'beds/*'), path(tsvs, stageAs: 'tsvs/*')

    output:
    tuple val(meta), path("orf_catalogue.bed12")   , emit: catalogue_bed
    tuple val(meta), path("orf_catalogue.tsv")     , emit: catalogue_tsv
    tuple val(meta), path("orf_to_gene.tsv")       , emit: orf_to_gene
    tuple val(meta), path("orf_catalogue.mqc.tsv") , emit: mqc
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/^Python //"'), topic: versions
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version 2>&1 | sed "s/^bedtools //"'), topic: versions
    tuple val("${task.process}"), val('cd-hit'), eval('cd-hit -h 2>&1 | head -n1 | sed "s/.*CD-HIT version //; s/ .*//"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p out
    orf_catalogue_merge.py \\
        --bed ${bed12s.collect{ it.toString() }.join(' ')} \\
        --tsv ${tsvs.collect{ it.toString() }.join(' ')} \\
        --out-dir out \\
        ${args}
    mv out/orf_catalogue.bed12 out/orf_catalogue.tsv out/orf_to_gene.tsv out/orf_catalogue.mqc.tsv .
    """

    stub:
    """
    touch orf_catalogue.bed12
    touch orf_catalogue.tsv
    touch orf_to_gene.tsv
    touch orf_catalogue.mqc.tsv
    """
}
