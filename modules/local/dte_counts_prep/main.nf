process DTE_COUNTS_PREP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1019bd22c111267bcb670fdb128829776f0ca6adfa7b0e2d126f91577d08e3/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:75514f9f977be607' }"

    input:
    tuple val(meta), path(primary_counts), path(secondary_counts), path(mapping)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: counts
    path "versions.yml"                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''
    template 'dte_counts_prep.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    python - <<END
import platform
import yaml

with open("versions.yml", "w") as fh:
    yaml.safe_dump(
        {
            "${task.process}": {
                "python": platform.python_version(),
                "pandas": "stub",
            }
        },
        fh,
        default_flow_style=False,
        sort_keys=False,
    )
END
    """
}
