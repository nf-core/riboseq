#!/usr/bin/env python3
"""Merge per-sample ORF P-site count TSVs into a single ORF x sample matrix.

Input files: any number of per-sample TSVs with three columns:
  sample_id  orf_id  count

The sample_id is identical on every row of a given file (it is prepended
upstream by QUANTIFY_INFRAME_PSITE_ORF). Files are passed as positional
arguments; the column order in the output matrix is the lexicographic sort
of the sample ids encountered.

Output: a single TSV with header `orf_id<tab>sample1<tab>sample2...` and
one row per ORF observed in the catalogue input (`--orf-list`). ORFs with
no counts in a given sample are filled with 0.

The `--orf-list` argument is a BED-like file whose 4th column is the ORF
id; this guarantees the matrix contains every ORF in the catalogue, not
just the subset with non-zero counts in at least one sample.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("counts", nargs="+", help="Per-sample count TSVs (sample_id, orf_id, count).")
    p.add_argument("--orf-list", required=True, help="BED/TSV whose 4th column is the canonical ORF id list.")
    p.add_argument("-o", "--output", default="orf_psite_counts.tsv", help="Output matrix path.")
    return p.parse_args(argv)


def read_orf_ids(path: Path) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            orf_id = parts[3]
            if orf_id in seen:
                continue
            seen.add(orf_id)
            ordered.append(orf_id)
    return ordered


def read_counts(path: Path) -> tuple[str | None, dict[str, int]]:
    counts: dict[str, int] = {}
    sample_id: str | None = None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            s, orf, c = parts[0], parts[1], parts[2]
            if sample_id is None:
                sample_id = s
            elif sample_id != s:
                # Mixed sample ids in one file is a producer bug, not user input.
                raise ValueError(f"{path}: encountered sample id '{s}' after '{sample_id}'")
            try:
                counts[orf] = counts.get(orf, 0) + int(c)
            except ValueError:
                counts[orf] = counts.get(orf, 0) + int(float(c))
    return sample_id, counts


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    orf_ids = read_orf_ids(Path(args.orf_list))
    per_sample: dict[str, dict[str, int]] = {}
    for path in args.counts:
        sample, counts = read_counts(Path(path))
        if sample is None:
            # Empty file - represent the sample with zeros if we can recover
            # its id from the filename, otherwise skip. Filename convention
            # used upstream is "<sample>.orf_inframe_psite_counts.tsv".
            stem = Path(path).name.split(".")[0]
            sample = stem
            counts = {}
        per_sample[sample] = counts

    samples = sorted(per_sample.keys())
    with open(args.output, "w") as out:
        out.write("orf_id\t" + "\t".join(samples) + "\n")
        for orf in orf_ids:
            row = [orf]
            for s in samples:
                row.append(str(per_sample[s].get(orf, 0)))
            out.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
