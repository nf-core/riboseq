#!/usr/bin/awk -f
#
# Build a full-reference + novel transcriptome annotation: every row of the
# full reference GTF (first file), then the rows of the hybrid GTF (second
# file) whose transcript the reference does not carry. Rows with no
# transcript_id are compared on gene_id.
#
# Consumers are the ORF catalogue and the ORF-level DTE RNA denominator, both
# of which need a superset of every annotation the ORF callers were given.
#
BEGIN { FS = OFS = "\t" }

FNR == NR {
    if ($0 !~ /^#/) {
        if (match($9, /transcript_id "([^"]+)"/, m)) ref_tx[m[1]]   = 1
        if (match($9, /gene_id "([^"]+)"/, m))       ref_gene[m[1]] = 1
    }
    print
    next
}

/^#/ { next }

match($9, /transcript_id "([^"]+)"/, m) { if (!(m[1] in ref_tx)) print; next }

match($9, /gene_id "([^"]+)"/, m) && !(m[1] in ref_gene)
