#!/usr/bin/awk -f
#
# Build a full-reference + novel transcriptome annotation for the ORF-level
# DTE denominator. Emits every row of the full reference GTF (first file),
# then the rows of the hybrid GTF (second file) whose gene_id is absent from
# the reference - i.e. the novel intergenic genes. Known-gene rows in the
# hybrid are already covered by the reference and are not duplicated.
#
BEGIN { FS = OFS = "\t" }

FNR == NR {
    if ($0 !~ /^#/ && match($9, /gene_id "([^"]+)"/, m)) ref_gene[m[1]] = 1
    print
    next
}

/^#/ { next }

match($9, /gene_id "([^"]+)"/, m) && !(m[1] in ref_gene)
