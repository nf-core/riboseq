#!/usr/bin/env python3
"""Aggregate per-ORF P-site counts to gene level using ONLY canonical CDS ORFs.

Produces a long-format `sample<TAB>gene_id<TAB>count` table matching the
shape of the per-sample `gene_inframe_psite_counts.tsv` files emitted by
`QUANTIFY_INFRAME_PSITE_PLASTID`, so the downstream
`REPLACE_RIBOSEQ_COUNTS_IN_MATRIX` awk substitution consumes it unchanged.
uORF / dORF / novel_u / smORF / other ORFs are excluded so the
translational numerator stays clean.

Inputs (positional, via the Nextflow template):
  ${orf_count_matrix}   ORF x sample matrix from issue #166.
                        First column `orf_id`, remaining columns are
                        sample ids, integer counts.
  ${orf_to_gene_tsv}    ORF id -> gene id mapping from issue #167.
                        Three columns: orf_id, gene_id, transcript_id.
                        First gene_id wins per ORF (deterministic by
                        file order) so an ORF's counts are never
                        double-counted at cohort level.
  ${orf_catalogue_tsv}  ORF catalogue sidecar TSV (carries
                        `orf_class`). Used to filter on `orf_class`.

`ext.args` (parsed below via argparse):
  --orf-class STR  ORF class to include in the gene sum.
                   Default: `canonical_cds`.

Output:
  `${prefix}.gene_cds_psite_counts.tsv`
"""
import argparse
import csv
import platform
import sys
from collections import defaultdict

import yaml


def load_canonical_orfs(catalogue_path, target_class):
    keep = set()
    with open(catalogue_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None or "orf_id" not in reader.fieldnames:
            raise SystemExit(
                f"{catalogue_path}: missing required `orf_id` header"
            )
        if "orf_class" not in reader.fieldnames:
            raise SystemExit(
                f"{catalogue_path}: missing required `orf_class` header"
            )
        for row in reader:
            if row.get("orf_class") == target_class:
                keep.add(row["orf_id"])
    return keep


def load_orf_to_gene(o2g_path):
    """First gene_id wins per ORF. Stable because the catalogue emits
    rows in cluster-iteration order, which is deterministic across runs."""
    mapping = {}
    with open(o2g_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if (
            reader.fieldnames is None
            or "orf_id" not in reader.fieldnames
            or "gene_id" not in reader.fieldnames
        ):
            raise SystemExit(
                f"{o2g_path}: expected `orf_id` and `gene_id` headers"
            )
        for row in reader:
            orf = row["orf_id"]
            gene = row.get("gene_id", "")
            if not orf or not gene:
                continue
            if orf in mapping:
                continue
            mapping[orf] = gene
    return mapping


def write_versions():
    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {
                "${task.process}": {
                    "python": platform.python_version(),
                }
            },
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orf-class",
        default="canonical_cds",
        help="ORF class to include in the gene sum. Default canonical_cds.",
    )
    raw_args = "${args}".split() if "${args}".strip() else []
    parsed_args = parser.parse_args(raw_args)

    keep_orfs = load_canonical_orfs("${orf_catalogue_tsv}", parsed_args.orf_class)
    orf_to_gene = load_orf_to_gene("${orf_to_gene_tsv}")

    with open("${orf_count_matrix}", newline="") as fh:
        reader = csv.reader(fh, delimiter="\\t")
        header = next(reader, None)
        if not header or header[0] != "orf_id":
            raise SystemExit(
                "${orf_count_matrix}: expected `orf_id` first-column header"
            )
        samples = header[1:]
        gene_sums = defaultdict(lambda: [0] * len(samples))
        for row in reader:
            if not row or not row[0]:
                continue
            orf_id = row[0]
            if orf_id not in keep_orfs:
                continue
            gene_id = orf_to_gene.get(orf_id)
            if not gene_id:
                continue
            for i, val in enumerate(row[1:]):
                try:
                    gene_sums[gene_id][i] += int(val)
                except ValueError:
                    gene_sums[gene_id][i] += int(float(val))

    with open("${prefix}.gene_cds_psite_counts.tsv", "w") as out:
        for gene_id in sorted(gene_sums.keys()):
            counts = gene_sums[gene_id]
            for i, c in enumerate(counts):
                out.write(f"{samples[i]}\\t{gene_id}\\t{c}\\n")

    write_versions()
    return 0


sys.exit(main())
