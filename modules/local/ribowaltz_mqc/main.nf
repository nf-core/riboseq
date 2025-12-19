process RIBOWALTZ_MQC {
    tag "ribowaltz_mqc"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    path psite_region_files
    path frames_files

    output:
    path "ribowaltz_psite_regions_mqc.tsv", emit: psite_regions_mqc, optional: true
    path "ribowaltz_frames_mqc.tsv"       , emit: frames_mqc      , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ribowaltz_to_mqc.py psite_region ${psite_region_files} > ribowaltz_psite_regions_mqc.tsv
    ribowaltz_to_mqc.py frames ${frames_files} > ribowaltz_frames_mqc.tsv
    """
}
