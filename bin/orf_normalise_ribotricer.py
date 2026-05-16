#!/usr/bin/env python3
"""Normalise a Ribotricer `detect-orfs` per-sample output to unified BED12 +
sidecar TSV.

Ribotricer writes `*_translating_ORFs.tsv` with columns:

    ORF_ID, ORF_type, status, chrom, strand, start_codon, gene_id, gene_name,
    gene_type, transcript_id, transcript_type, length, profile, coordinate,
    valid_codons, valid_codons_ratio, read_count, count, count_ratio,
    phase_score, codon_idx, used_lengths

The `coordinate` column lists exon-segment ranges (1-based inclusive),
comma-separated, e.g. '12345-12500,12750-12900'. We use it directly to
build BED12 blocks; this means the GTF lookup is only needed for sanity
filtering.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_normalise_common import (
    classify_ribotricer,
    emit_bed12,
    emit_tsv_row,
    load_transcripts,
    reclassify_smorf,
    write_outputs,
)


def parse_coordinate(s: str) -> list[tuple[int, int]]:
    """Parse Ribotricer 'coordinate' field. Returns list of genomic blocks
    (0-based half-open), sorted ascending by start."""
    if not s:
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok or "-" not in tok:
            continue
        a, b = tok.split("-", 1)
        try:
            a_i, b_i = int(a), int(b)
        except ValueError:
            continue
        if b_i < a_i:
            a_i, b_i = b_i, a_i
        out.append((a_i - 1, b_i))
    out.sort()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    args = ap.parse_args()

    if not args.input.exists() or args.input.stat().st_size == 0:
        write_outputs(args.out_bed, args.out_tsv, [], [])
        return 0

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
            status = (row.get("status") or "").lower()
            if status and status != "translating":
                continue
            orf_id_raw = row.get("ORF_ID") or ""
            tid = row.get("transcript_id") or ""
            gid = row.get("gene_id") or ""
            chrom = row.get("chrom") or ""
            strand = row.get("strand") or "+"
            try:
                length_nt = int(row.get("length") or 0)
            except ValueError:
                length_nt = 0
            # ORF length is nt of full ORF including stop; aa_length = (length - 3)/3
            aa_len = max(0, (length_nt - 3) // 3) if length_nt > 0 else 0
            orf_type = row.get("ORF_type") or ""
            orf_class = classify_ribotricer(orf_type)
            orf_class = reclassify_smorf(orf_class, aa_len)

            # Ribotricer's phase_score - higher = better. Map 0..1 to 0..1000.
            score_raw = row.get("phase_score") or "0"
            try:
                pscore = float(score_raw)
                score_int = max(0, min(1000, int(round(pscore * 1000))))
            except (TypeError, ValueError):
                pscore = float("nan")
                score_int = 0

            blocks = parse_coordinate(row.get("coordinate") or "")
            if not blocks:
                continue

            orf_id = f"ribotricer|{orf_id_raw}"
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
                    caller="ribotricer",
                    sample_id=args.sample_id,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    gene_id=gid,
                    transcript_id=tid,
                    orf_class=orf_class,
                    aa_length=aa_len,
                    score="",  # ribotricer score excluded from ranking (#163)
                )
            )

    write_outputs(args.out_bed, args.out_tsv, bed_lines, tsv_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
