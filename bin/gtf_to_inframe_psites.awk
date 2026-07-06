#!/usr/bin/awk -f
# Convert GTF CDS segments to a BED file of in-frame p-site positions.
#
# For each CDS entry in the GTF file, this script calculates the positions
# of codon starts (i.e., potential in-frame p-sites) based on the phase
# annotation. The output is a deduplicated BED6 file with columns:
#   chrom  start  end  feature_id  0  strand
#
# Usage:
#   awk -v FEATURE=gene -f gtf_to_inframe_psites.awk annotation.gtf > inframe_psites.bed
#
# The FEATURE variable controls which GTF attribute is used to group counts
# (e.g., "gene" extracts gene_id, "transcript" extracts transcript_id).

BEGIN { FS = OFS = "\t" }

$3 == "CDS" && $9 ~ (FEATURE "_id \"") {

    # parse GTF fields
    chrom = $1; start = $4; end = $5; strand = $7; phase = $8
    feature = $9; sub(".*" FEATURE "_id \"", "", feature); sub(/".*/, "", feature)

    # print deduplicated coordinates of codon starts
    frame = strand == "+" ? phase : (3 - phase + end - start) % 3
    for (pos = start + frame; pos <= end; pos += 3) {
        $0 = chrom OFS pos-1 OFS pos OFS feature OFS 0 OFS strand
        if (!seen[$0]++) print
    }
}
