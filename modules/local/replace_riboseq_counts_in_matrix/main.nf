// Replace Ribo-seq counts in a count matrix with in-frame p-site counts.
// The p-site counts are expected as a three-column TSV (sample, feature, count),
// e.g., as produced by QUANTIFY_INFRAME_PSITE_PLASTID. For every Ribo-seq sample found
// in the p-site counts, the corresponding column in the count matrix is updated.

process REPLACE_RIBOSEQ_COUNTS_IN_MATRIX {
    tag "$counts"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(psite_counts)
    tuple val(meta2), path(counts)
    val feature

    output:
    tuple val(meta2), path("*_counts.tsv"), emit: counts
    path "versions.yml", emit: versions

    script:
    """
    awk -F'\\t' -v OFS='\\t' '

        # receive the new counts from stdin
        FILENAME == "/dev/stdin" {
            # cache sample names and counts
            samples[\$1]
            count[\$1, \$2] = \$3
        }

        # the old count matrix comes after stdin
        FILENAME != "/dev/stdin" {
            if (FNR == 1) {
                # the column headers tell us the sample names
                for (i = 3; i <= NF; i++)
                    header[i] = \$i
            } else {
                # replace values whenever we find a column of a Ribo-seq sample
                for (i = 3; i <= NF; i++)
                    if (header[i] in samples)
                        \$i = count[header[i], \$1] + 0
            }
            print
        }
    ' /dev/stdin ${counts} < ${psite_counts} > "${feature}_counts.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -f3 -d' ' | cut -f1 -d,)
    END_VERSIONS
    """
}
