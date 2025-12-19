#!/usr/bin/env python3
"""
Transform ribowaltz TSV output to MultiQC-compatible wide format.

Usage:
    ribowaltz_to_mqc.py psite_region <input_files...> > output.tsv
    ribowaltz_to_mqc.py frames <input_files...> > output.tsv
"""

import sys
import csv
from collections import defaultdict


def transform_psite_region(files):
    """Transform psite_region.tsv files to wide format for MultiQC bargraph."""
    data = defaultdict(dict)
    regions = ["5' UTR", "CDS", "3' UTR"]

    for filepath in files:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample = row['sample']
                region = row['region']
                # Skip the 'RNAs' row which is the reference distribution
                if sample == 'RNAs':
                    continue
                # Map region names to cleaner format
                region_map = {"5' UTR": "5' UTR", "5utr": "5' UTR",
                              "CDS": "CDS", "cds": "CDS",
                              "3' UTR": "3' UTR", "3utr": "3' UTR"}
                clean_region = region_map.get(region, region)
                if clean_region in regions:
                    data[sample][clean_region] = float(row['scaled_count'])

    # Output wide format
    print("Sample\t" + "\t".join(regions))
    for sample in sorted(data.keys()):
        values = [str(round(data[sample].get(r, 0), 1)) for r in regions]
        print(f"{sample}\t" + "\t".join(values))


def transform_frames(files):
    """Transform frames.tsv files to MultiQC custom content.

    Structure matches ribowaltz frames.pdf plot:
    - Groups: regions (5' UTR, CDS, 3' UTR) - using sample_groups
    - Within each group: one bar per sample
    - Each bar stacked: Frame 0, Frame 1, Frame 2 segments

    Key QC insight: CDS should show Frame 0 enrichment, UTRs should not.
    """
    data = defaultdict(lambda: defaultdict(dict))
    frames = ["Frame 0", "Frame 1", "Frame 2"]
    regions = ["5utr", "cds", "3utr"]
    region_labels = {"5utr": "5' UTR", "cds": "CDS", "3utr": "3' UTR"}

    region_map = {
        "5' UTR": "5utr", "5utr": "5utr",
        "CDS": "cds", "cds": "cds",
        "3' UTR": "3utr", "3utr": "3utr"
    }

    for filepath in files:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample = row['sample']
                region = region_map.get(row['region'], row['region'])
                frame = f"Frame {row['frame']}"
                if region in regions and frame in frames:
                    data[sample][region][frame] = float(row['scaled_count'])

    samples = sorted(data.keys())

    # Output YAML header
    print("# id: 'ribowaltz_frames'")
    print("# section_name: 'riboWaltz: Reading Frame Distribution'")
    print("# description: 'Distribution of P-sites across reading frames for each transcript region. Good Ribo-seq data shows Frame 0 enrichment (>50%) in the CDS but not in UTRs.'")
    print("# plot_type: 'bargraph'")
    print("# pconfig:")
    print("#     id: 'ribowaltz_frames_plot'")
    print("#     title: 'riboWaltz: Reading Frame Distribution'")
    print("#     ylab: '% of P-sites'")
    print("#     cpswitch: false")

    # Build sample_groups: group by region, with samples within each region group
    # Structure: {"5' UTR": [["sample1_5' UTR", "sample1"], ...], "CDS": [...], ...}
    if len(samples) > 0:
        print("#     sample_groups:")
        for region in regions:
            label = region_labels[region]
            print(f'#         "{label}":')
            for sample in samples:
                # [sample_key, display_label] - display_label is the sample name
                print(f'#             - ["{sample}_{label}", "{sample}"]')

    # Output data: each row is {sample}_{region} with frame values
    print("Sample\t" + "\t".join(frames))
    for sample in samples:
        for region in regions:
            label = region_labels[region]
            row_id = f"{sample}_{label}"
            values = [str(round(data[sample][region].get(f, 0), 1)) for f in frames]
            print(f"{row_id}\t" + "\t".join(values))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    mode = sys.argv[1]
    files = sys.argv[2:]

    if mode == "psite_region":
        transform_psite_region(files)
    elif mode == "frames":
        transform_frames(files)
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        print(__doc__, file=sys.stderr)
        sys.exit(1)
