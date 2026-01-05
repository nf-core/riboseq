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

# Common constants
PARENT_ID = "ribowaltz"
PARENT_NAME = "riboWaltz"
PARENT_DESC = 'Quality control metrics from <a href="https://github.com/LabTranslationalArchitectomics/riboWaltz">riboWaltz</a> for assessing Ribo-seq data quality.'

REGIONS = ["5' UTR", "CDS", "3' UTR"]
REGION_MAP = {"5' UTR": "5' UTR", "5utr": "5' UTR", "CDS": "CDS", "cds": "CDS", "3' UTR": "3' UTR", "3utr": "3' UTR"}


def read_tsv_files(files):
    """Yield rows from multiple TSV files."""
    for filepath in files:
        with open(filepath, 'r') as f:
            yield from csv.DictReader(f, delimiter='\t')


def print_bargraph_header(plot_id, section_name, description, ylab="% of P-sites"):
    """Print common YAML header for bargraph plots."""
    print(f"# id: '{plot_id}'")
    print(f"# section_name: '{section_name}'")
    print(f"# description: '{description}'")
    print(f"# parent_id: '{PARENT_ID}'")
    print(f"# parent_name: '{PARENT_NAME}'")
    print(f"# parent_description: '{PARENT_DESC}'")
    print("# plot_type: 'bargraph'")
    print("# pconfig:")
    print(f"#     id: '{plot_id}_plot'")
    print(f"#     title: '{PARENT_NAME}: {section_name}'")
    print(f"#     ylab: '{ylab}'")
    print("#     cpswitch: false")


def transform_psite_region(files):
    """Transform psite_region.tsv files to wide format for MultiQC bargraph."""
    data = defaultdict(dict)
    rna_ref = {}

    for row in read_tsv_files(files):
        sample = row['sample']
        region = REGION_MAP.get(row['region'], row['region'])

        if sample == 'RNAs':
            if region in REGIONS and region not in rna_ref:
                rna_ref[region] = float(row['scaled_count'])
        elif region in REGIONS:
            data[sample][region] = float(row['scaled_count'])

    print_bargraph_header(
        "ribowaltz_1_psite_regions",
        "P-site Region Distribution",
        "Distribution of P-sites across transcript regions. Good Ribo-seq data shows strong CDS enrichment (>70%). The RNA-seq reference shows expected distribution from uniform transcript coverage."
    )

    print("Sample\t" + "\t".join(REGIONS))
    if rna_ref:
        values = [str(round(rna_ref.get(r, 0), 1)) for r in REGIONS]
        print(f"<b>RNA-seq reference</b>\t" + "\t".join(values))
    for sample in sorted(data.keys()):
        values = [str(round(data[sample].get(r, 0), 1)) for r in REGIONS]
        print(f"{sample}\t" + "\t".join(values))


def transform_frames(files):
    """Transform frames.tsv files to MultiQC custom content."""
    data = defaultdict(lambda: defaultdict(dict))
    frames = ["Frame 0", "Frame 1", "Frame 2"]
    region_keys = ["5utr", "cds", "3utr"]
    region_labels = {"5utr": "5' UTR", "cds": "CDS", "3utr": "3' UTR"}
    to_key = {"5' UTR": "5utr", "5utr": "5utr", "CDS": "cds", "cds": "cds", "3' UTR": "3utr", "3utr": "3utr"}

    for row in read_tsv_files(files):
        sample = row['sample']
        region = to_key.get(row['region'], row['region'])
        frame = f"Frame {row['frame']}"
        if region in region_keys and frame in frames:
            data[sample][region][frame] = float(row['scaled_count'])

    samples = sorted(data.keys())

    print_bargraph_header(
        "ribowaltz_2_frames",
        "Reading Frame Distribution",
        "Distribution of P-sites across reading frames for each transcript region. Good Ribo-seq data shows Frame 0 enrichment (>50%) in the CDS but not in UTRs."
    )

    if samples:
        print("#     sample_groups:")
        for region in region_keys:
            label = region_labels[region]
            print(f'#         "{label}":')
            for sample in samples:
                print(f'#             - ["{sample}_{label}", "{sample}"]')

    print("Sample\t" + "\t".join(frames))
    for sample in samples:
        for region in region_keys:
            label = region_labels[region]
            values = [str(round(data[sample][region].get(f, 0), 1)) for f in frames]
            print(f"{sample}_{label}\t" + "\t".join(values))


def transform_metaprofile(files, region_filter):
    """Transform metaprofile_psite.tsv files to MultiQC linegraph JSON format."""
    data = defaultdict(list)
    target_region = f"Distance from {region_filter} (nt)"

    for row in read_tsv_files(files):
        if row['region'] == target_region:
            data[row['sample']].append([int(float(row['x'])), round(float(row['y']), 2)])

    for sample in data:
        data[sample].sort(key=lambda p: p[0])

    section_name = f"Metaprofile ({region_filter.title()} Codon)"
    plot_id = f"ribowaltz_{'3' if region_filter == 'start' else '4'}_metaprofile_{region_filter}"

    output = {
        "id": plot_id,
        "section_name": section_name,
        "description": f"P-site frequency around the {region_filter} codon. Good Ribo-seq data shows trinucleotide periodicity with peaks at frame 0 positions.",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "parent_description": PARENT_DESC,
        "plot_type": "linegraph",
        "pconfig": {
            "id": f"{plot_id}_plot",
            "title": f"{PARENT_NAME}: {section_name}",
            "xlab": f"Distance from {region_filter} codon (nt)",
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

    mode, files = sys.argv[1], sys.argv[2:]
    modes = {
        "psite_region": lambda: transform_psite_region(files),
        "frames": lambda: transform_frames(files),
        "metaprofile_start": lambda: transform_metaprofile(files, "start"),
        "metaprofile_stop": lambda: transform_metaprofile(files, "stop"),
    }

    if mode in modes:
        modes[mode]()
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        print(__doc__, file=sys.stderr)
        sys.exit(1)
