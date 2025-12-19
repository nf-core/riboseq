#!/usr/bin/env python3
"""
Transform ribowaltz TSV output to MultiQC-compatible format.

Usage:
    ribowaltz_to_mqc.py psite_region <input_files...> > output.tsv
    ribowaltz_to_mqc.py frames <input_files...> > output.tsv
    ribowaltz_to_mqc.py metaprofile_start <input_files...> > output.json
    ribowaltz_to_mqc.py metaprofile_stop <input_files...> > output.json
"""

import sys
import csv
import json
from collections import defaultdict


def transform_psite_region(files):
    """Transform psite_region.tsv files to wide format for MultiQC bargraph."""
    data = defaultdict(dict)
    rna_ref = {}  # Global RNA-seq reference from annotation
    regions = ["5' UTR", "CDS", "3' UTR"]

    region_map = {"5' UTR": "5' UTR", "5utr": "5' UTR",
                  "CDS": "CDS", "cds": "CDS",
                  "3' UTR": "3' UTR", "3utr": "3' UTR"}

    for filepath in files:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample = row['sample']
                region = row['region']
                clean_region = region_map.get(region, region)

                # Capture the 'RNAs' reference (same for all samples using same annotation)
                if sample == 'RNAs':
                    if clean_region in regions and clean_region not in rna_ref:
                        rna_ref[clean_region] = float(row['scaled_count'])
                    continue

                if clean_region in regions:
                    data[sample][clean_region] = float(row['scaled_count'])

    samples = sorted(data.keys())

    # Output YAML header
    print("# id: 'ribowaltz_1_psite_regions'")
    print("# section_name: 'P-site Region Distribution'")
    print("# description: 'Distribution of P-sites across transcript regions. Good Ribo-seq data shows strong CDS enrichment (>70%). The RNA-seq reference shows expected distribution from uniform transcript coverage.'")
    print("# parent_id: 'ribowaltz'")
    print("# parent_name: 'riboWaltz'")
    print("# parent_description: 'Quality control metrics from <a href=\"https://github.com/LabTranslationalArchitectomics/riboWaltz\">riboWaltz</a> for assessing Ribo-seq data quality.'")
    print("# plot_type: 'bargraph'")
    print("# pconfig:")
    print("#     id: 'ribowaltz_1_psite_regions_plot'")
    print("#     title: 'riboWaltz: P-site Region Distribution'")
    print("#     ylab: '% of P-sites'")
    print("#     cpswitch: false")

    # Output data
    print("Sample\t" + "\t".join(regions))

    # Output RNA-seq reference first if available
    if rna_ref:
        values = [str(round(rna_ref.get(r, 0), 1)) for r in regions]
        print(f"<b>RNA-seq reference</b>\t" + "\t".join(values))

    for sample in samples:
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
    print("# id: 'ribowaltz_2_frames'")
    print("# section_name: 'Reading Frame Distribution'")
    print("# description: 'Distribution of P-sites across reading frames for each transcript region. Good Ribo-seq data shows Frame 0 enrichment (>50%) in the CDS but not in UTRs.'")
    print("# parent_id: 'ribowaltz'")
    print("# parent_name: 'riboWaltz'")
    print("# parent_description: 'Quality control metrics from <a href=\"https://github.com/LabTranslationalArchitectomics/riboWaltz\">riboWaltz</a> for assessing Ribo-seq data quality.'")
    print("# plot_type: 'bargraph'")
    print("# pconfig:")
    print("#     id: 'ribowaltz_2_frames_plot'")
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


def transform_metaprofile(files, region_filter):
    """Transform metaprofile_psite.tsv files to MultiQC linegraph JSON format.

    Args:
        files: List of metaprofile TSV files
        region_filter: "start" or "stop" to filter by region type

    The input TSV has columns: sample, region, x, y
    Region values are "Distance from start (nt)" or "Distance from stop (nt)"
    """
    data = defaultdict(list)

    region_map = {
        "start": "Distance from start (nt)",
        "stop": "Distance from stop (nt)"
    }
    target_region = region_map[region_filter]

    for filepath in files:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample = row['sample']
                region = row['region']
                if region == target_region:
                    x = int(float(row['x']))
                    y = round(float(row['y']), 2)
                    data[sample].append([x, y])

    # Sort each sample's data by x value
    for sample in data:
        data[sample].sort(key=lambda point: point[0])

    # Build JSON output
    if region_filter == "start":
        section_name = "Metaprofile (Start Codon)"
        description = "P-site frequency around the start codon. Good Ribo-seq data shows trinucleotide periodicity with peaks at frame 0 positions."
        plot_id = "ribowaltz_3_metaprofile_start"
        xlab = "Distance from start codon (nt)"
    else:
        section_name = "Metaprofile (Stop Codon)"
        description = "P-site frequency around the stop codon. Good Ribo-seq data shows trinucleotide periodicity with peaks at frame 0 positions."
        plot_id = "ribowaltz_4_metaprofile_stop"
        xlab = "Distance from stop codon (nt)"

    output = {
        "id": plot_id,
        "section_name": section_name,
        "description": description,
        "parent_id": "ribowaltz",
        "parent_name": "riboWaltz",
        "parent_description": "Quality control metrics from <a href=\"https://github.com/LabTranslationalArchitectomics/riboWaltz\">riboWaltz</a> for assessing Ribo-seq data quality.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": f"{plot_id}_plot",
            "title": f"riboWaltz: {section_name}",
            "xlab": xlab,
            "ylab": "P-site frequency",
            "x_decimals": False
        },
        "data": dict(sorted(data.items()))
    }

    print(json.dumps(output, indent=2))


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
    elif mode == "metaprofile_start":
        transform_metaprofile(files, "start")
    elif mode == "metaprofile_stop":
        transform_metaprofile(files, "stop")
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        print(__doc__, file=sys.stderr)
        sys.exit(1)
