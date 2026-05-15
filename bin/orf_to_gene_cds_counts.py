#!/usr/bin/env python3
"""Aggregate per-ORF P-site counts to gene level using ONLY canonical CDS ORFs.

This is the Tier 1 input-replacement helper for issue #168. It produces a
gene x sample count matrix that can be substituted into the existing
gene-level translational-efficiency input matrix in place of the
plastid-derived gene-level CDS p-site counts.

Inputs:
  --orf-counts orf_psite_counts.tsv
    ORF x sample matrix from issue #166. First column `orf_id`, remaining
    columns are sample ids, integer counts.
  --orf-to-gene orf_to_gene.tsv
    Mapping from issue #167. Three columns: orf_id, gene_id, transcript_id.
    One ORF can map to multiple (gene_id, transcript_id) rows when callers
    placed it on different annotated isoforms; this script collapses by
    gene_id, picking the first gene_id observed per ORF (deterministic via
    file order) so an ORF's counts are never double-counted at the
    cohort level.
  --orf-catalogue orf_catalogue.tsv
    Per-ORF classification table from issue #167. Used to filter to
    `orf_class == "canonical_cds"`; uORFs, dORFs, novel_u, smORFs and
    "other" are excluded from the gene sum so that the translational
    numerator stays clean (see issue #168 spec, Tier 1 rationale).

Output:
  Long-format TSV with three columns (no header), one row per
  (sample, gene_id) with the integer sum of canonical-CDS ORF P-site
  counts:
    <sample><TAB><gene_id><TAB><count>
  This matches the shape of the per-sample
  `gene_inframe_psite_counts.tsv` files produced by
  QUANTIFY_INFRAME_PSITE_PLASTID, so the same downstream
  REPLACE_RIBOSEQ_COUNTS_IN_MATRIX awk script consumes it without
  modification.

Notes:
  - Genes with no canonical-CDS ORF in the catalogue are NOT emitted.
    This is consistent with the existing gene-level plastid input, which
    also only carries rows for genes with at least one in-frame P-site.
  - Sample columns are emitted in lexicographic order matching the input
    matrix header.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--orf-counts", required=True, type=Path)
    p.add_argument("--orf-to-gene", required=True, type=Path)
    p.add_argument("--orf-catalogue", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument(
        "--orf-class",
        default="canonical_cds",
        help="ORF class to include in the gene sum. Default canonical_cds.",
    )
    return p.parse_args(argv)


def load_canonical_orfs(catalogue_path: Path, target_class: str) -> set[str]:
    keep: set[str] = set()
    with open(catalogue_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
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


def load_orf_to_gene(o2g_path: Path) -> dict[str, str]:
    """First gene_id wins per ORF. Stable because the catalogue emits
    rows in cluster-iteration order, which is deterministic across runs."""
    mapping: dict[str, str] = {}
    with open(o2g_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "orf_id" not in reader.fieldnames or "gene_id" not in reader.fieldnames:
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


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    keep_orfs = load_canonical_orfs(args.orf_catalogue, args.orf_class)
    orf_to_gene = load_orf_to_gene(args.orf_to_gene)

    with open(args.orf_counts, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "orf_id":
            raise SystemExit(
                f"{args.orf_counts}: expected `orf_id` first-column header"
            )
        samples = header[1:]
        # gene -> [count per sample]
        gene_sums: dict[str, list[int]] = defaultdict(lambda: [0] * len(samples))
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

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        for gene_id in sorted(gene_sums.keys()):
            counts = gene_sums[gene_id]
            for i, c in enumerate(counts):
                out.write(f"{samples[i]}\t{gene_id}\t{c}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
