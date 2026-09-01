#!/usr/bin/awk -f
# Filter a gene-keyed counts TSV down to the gene_id set found in a
# canonical GTF. The header row is preserved; data rows are emitted only
# when column 1 (gene_id) appears as a `gene` or `transcript` feature in
# the GTF.
#
# Usage (GTF first, counts second):
#   awk -f filter_counts_canonical.awk canonical.gtf counts.tsv > canonical_filtered.tsv

BEGIN { FS = OFS = "\t" }

# First file (canonical GTF): collect gene_ids from gene / transcript rows.
FNR == NR {
    if ($0 ~ /^#/) { next }
    if ($3 != "gene" && $3 != "transcript") { next }
    if (match($9, /gene_id "([^"]+)"/, m)) keep[m[1]] = 1
    next
}

# Second file (counts TSV): preserve header, retain rows whose gene_id is kept.
FNR == 1 { print; next }
($1 in keep) { print }
