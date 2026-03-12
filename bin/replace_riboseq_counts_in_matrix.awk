#!/usr/bin/awk -f
# Replace Ribo-seq counts in a count matrix with in-frame p-site counts.
#
# The input file whose name matches "inframe_psite_counts" contains the new
# (p-site) counts as three columns:
#   sample \t feature \t count
#
# The other input file is the existing count matrix (with a header row).
# For every sample found in the p-site counts, the corresponding column
# in the matrix is replaced.

BEGIN { FS = OFS = "\t" }

# P-site counts file
FILENAME ~ /inframe_psite_counts/ {
    samples[$1]
    count[$1, $2] = $3
}

# Count matrix
FILENAME !~ /inframe_psite_counts/ {
    if (FNR == 1) {
        # The column headers tell us the sample names
        for (i = 3; i <= NF; i++)
            header[i] = $i
    } else {
        # Replace values whenever we find a column of a Ribo-seq sample
        for (i = 3; i <= NF; i++)
            if (header[i] in samples)
                $i = count[header[i], $1] + 0
    }
    print
}
