// Count Ribo-seq P-sites that are in-frame for each ORF in the cohort
// catalogue. P-site density is supplied as per-sample forward / reverse
// bedGraph tracks (e.g. from PLASTID_MAKE_WIGGLE). In-frame codon-start
// positions are pre-computed once per cohort from the BED12 catalogue and
// passed in via `inframe_psites`.
//
// Output: one TSV per sample with three columns:
//   sample_id  orf_id  count
//
// This module is the ORF-level counterpart of QUANTIFY_INFRAME_PSITE_PLASTID
// (gene-level). The two run side-by-side: gene-level retains its existing
// role in the Salmon-matrix replacement path; ORF-level feeds the
// downstream DTE matrix.

process QUANTIFY_INFRAME_PSITE_ORF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools==2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(meta), path(forward_bedgraph), path(reverse_bedgraph)
    tuple val(meta2), path(inframe_psites)

    output:
    tuple val(meta), path("*.orf_inframe_psite_counts.tsv"), emit: counts
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # append strand to bedgraph rows and stream both files together
    (
        sed 's/\$/\\t0\\t+/' "$forward_bedgraph"
        sed 's/\$/\\t0\\t-/' "$reverse_bedgraph"
    ) |

    # intersect codon-start positions with stranded P-site density,
    # then aggregate counts per ORF id (column 4 of the inframe psites BED)
    bedtools intersect -wa -wb -a $inframe_psites -b /dev/stdin -s |
    sort -k4,4 |
    bedtools groupby -g 4 -c 10 -o sum |
    sed 's/^/${meta.id}\\t/' \\
    > "${meta.id}.orf_inframe_psite_counts.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | cut -f2 -dv)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.orf_inframe_psite_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | cut -f2 -dv)
    END_VERSIONS
    """
}
