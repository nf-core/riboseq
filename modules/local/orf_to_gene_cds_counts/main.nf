// Re-aggregate per-ORF P-site counts to gene level using ONLY canonical_cds
// ORFs. Output is long-format `sample<TAB>gene_id<TAB>count`, matching the
// shape of the per-sample plastid `gene_inframe_psite_counts.tsv` so the
// REPLACE_RIBOSEQ_COUNTS_IN_MATRIX awk substitution step consumes it without
// changes. uORF / dORF / novel_u / smORF / other ORFs are excluded from the
// sum to keep the translational numerator clean.

process ORF_TO_GENE_CDS_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1019bd22c111267bcb670fdb128829776f0ca6adfa7b0e2d126f91577d08e3/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:75514f9f977be607' }"

    input:
    tuple val(meta), path(orf_count_matrix), path(orf_to_gene_tsv), path(orf_catalogue_tsv)

    output:
    tuple val(meta), path("${prefix}.gene_cds_psite_counts.tsv"), emit: gene_counts
    path "versions.yml"                                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''
    template 'orf_to_gene_cds_counts.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_cds_psite_counts.tsv

    python - <<END
import platform
import yaml

with open("versions.yml", "w") as fh:
    yaml.safe_dump(
        {
            "${task.process}": {
                "python": platform.python_version(),
            }
        },
        fh,
        default_flow_style=False,
        sort_keys=False,
    )
END
    """
}
