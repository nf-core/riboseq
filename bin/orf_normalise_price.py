#!/usr/bin/env python3
"""Normalise PRICE's `${prefix}.orfs.tsv` to unified BED12 + sidecar TSV.

PRICE writes `${prefix}.orfs.tsv` with columns:

    Gene, Id, Location, Candidate Location, Codon, Type, Start, Range,
    p value, <Condition-name(s)>, Total

The `Location` field encodes the ORF block structure in Gedi notation:

    chrom:strand:start-end[|start-end]...

Coordinates are 0-based, start inclusive, end exclusive (matches BED).
"strand" is '+' or '-'. Multiple exon segments are pipe-separated.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_normalise_common import (
    classify_price,
    emit_bed12,
    emit_tsv_row,
    load_transcripts,
    reclassify_smorf,
    write_outputs,
)


_LOCATION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<strand>[+-]):(?P<blocks>.+)$")


def parse_location(loc: str) -> tuple[str, str, list[tuple[int, int]]]:
    """Parse a Gedi Location field into (chrom, strand, blocks).

    Blocks are 0-based half-open genomic intervals, sorted ascending by start.
    Returns ('', '+', []) when the field cannot be parsed.
    """
    if not loc:
        return "", "+", []
    m = _LOCATION_RE.match(loc.strip())
    if not m:
        return "", "+", []
    chrom = m.group("chrom")
    strand = m.group("strand")
    blocks: list[tuple[int, int]] = []
    for tok in m.group("blocks").split("|"):
        tok = tok.strip()
        if not tok or "-" not in tok:
            continue
        a, b = tok.split("-", 1)
        try:
            start = int(a)
            end = int(b)
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        blocks.append((start, end))
    blocks.sort()
    return chrom, strand, blocks


def parse_total(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_p_value(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def score_from_p(p: float | None) -> int:
    """Map a p value to a 0-1000 BED score. Lower p -> higher score."""
    if p is None or p <= 0:
        return 1000
    if p >= 1:
        return 0
    import math
    s = int(round(min(-math.log10(p) * 100, 1000)))
    return max(0, min(1000, s))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path,
                    help="PRICE ${prefix}.orfs.tsv")
    ap.add_argument("--gtf", required=True, type=Path,
                    help="Reference GTF (used only for the transcript lookup pass)")
    ap.add_argument("--sample-id", required=True,
                    help="Cohort label written into the sidecar TSV (PRICE is a "
                         "cohort-level caller).")
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    args = ap.parse_args()

    if not args.input.exists() or args.input.stat().st_size == 0:
        write_outputs(args.out_bed, args.out_tsv, [], [])
        return 0

    # PRICE doesn't emit gene_id in its TSV for some ORF classes (notably
    # intronic/orphan), but reports transcript_id-style ORF ids. We don't
    # strictly need the transcript lookup here, but load it anyway so that
    # the sidecar transcript_id column matches the rest of the catalogue.
    transcripts = load_transcripts(args.gtf)

    bed_lines: list[str] = []
    tsv_rows: list[str] = []
    seen_orf_ids: set[str] = set()

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            write_outputs(args.out_bed, args.out_tsv, [], [])
            return 0

        for row in reader:
            orf_id_raw = (row.get("Id") or "").strip()
            if not orf_id_raw:
                continue
            location = row.get("Location") or ""
            chrom, strand, blocks = parse_location(location)
            if not blocks:
                continue

            orf_type = (row.get("Type") or "").strip()
            orf_class = classify_price(orf_type)

            # PRICE's "Id" is structured "transcript_id__ORFType__N"; recover
            # transcript_id from it for the sidecar.
            tid = orf_id_raw.split("__", 1)[0]
            gid = (row.get("Gene") or "").strip()
            if not gid:
                t = transcripts.get(tid)
                if t is not None:
                    gid = t.gene_id

            length_nt = sum(b[1] - b[0] for b in blocks)
            # PRICE reports start-to-stop ranges that INCLUDE the stop codon.
            aa_len = max(0, (length_nt - 3) // 3) if length_nt > 0 else 0
            orf_class = reclassify_smorf(orf_class, aa_len)

            pval = parse_p_value(row.get("p value") or row.get("p_value") or "")
            score_int = score_from_p(pval)

            orf_id = f"price|{orf_id_raw}"
            if orf_id in seen_orf_ids:
                continue
            seen_orf_ids.add(orf_id)

            bed_lines.append(
                emit_bed12(
                    chrom=chrom,
                    blocks=blocks,
                    name=orf_id,
                    score=score_int,
                    strand=strand,
                )
            )
            tsv_rows.append(
                emit_tsv_row(
                    orf_id=orf_id,
                    caller="price",
                    sample_id=args.sample_id,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    gene_id=gid,
                    transcript_id=tid,
                    orf_class=orf_class,
                    aa_length=aa_len,
                    score=pval if pval is not None else "",
                )
            )

    write_outputs(args.out_bed, args.out_tsv, bed_lines, tsv_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
