#!/usr/bin/env python3
"""Normalise an Rp-Bp `predict-translated-orfs` per-sample BED12 output to
unified BED12 + sidecar TSV.

Rp-Bp emits `*.predicted-orfs.bed.gz`: a standard BED12 (12 columns)
followed by extra columns:

    chrom start end name score strand thickStart thickEnd itemRgb
    blockCount blockSizes blockStarts orf_num orf_len orf_type
    bayes_factor_mean ...

The `name` field encodes 'transcript_id_orf-num' (Rp-Bp 4.x). We rely on
the BED12 columns as-is and parse the extras to populate the sidecar.
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_normalise_common import (
    classify_rpbp,
    emit_bed12,
    emit_tsv_row,
    load_transcripts,
    reclassify_smorf,
    write_outputs,
)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


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

    # Ribosome-pruned context: load transcripts for transcript->gene lookups.
    transcripts = load_transcripts(args.gtf)

    bed_lines: list[str] = []
    tsv_rows: list[str] = []
    seen_orf_ids: set[str] = set()

    with open_text(args.input) as fh:
        # Determine columns from header if present, otherwise positional.
        first = True
        header: list[str] = []
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                if line.startswith("#") and first:
                    header = line.lstrip("#").strip().split("\t")
                    first = False
                continue
            first = False
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            name = cols[3]
            strand = cols[5] if cols[5] in ("+", "-") else "+"
            try:
                block_count = int(cols[9])
            except ValueError:
                continue
            block_sizes_raw = cols[10].rstrip(",").split(",")
            block_starts_raw = cols[11].rstrip(",").split(",")
            try:
                block_sizes = [int(x) for x in block_sizes_raw if x]
                block_starts = [int(x) for x in block_starts_raw if x]
            except ValueError:
                continue
            if not block_sizes or len(block_sizes) != block_count:
                continue
            blocks = [(start + bs, start + bs + sz) for bs, sz in zip(block_starts, block_sizes)]
            blocks.sort()

            # Extra columns (positional). Rp-Bp 4.x layout typically:
            #   12=orf_num, 13=orf_len, 14=orf_type, 15=bayes_factor_mean, ...
            orf_type = cols[14] if len(cols) > 14 else ""
            try:
                orf_len_nt = int(cols[13]) if len(cols) > 13 else 0
            except ValueError:
                orf_len_nt = 0
            try:
                bf_mean = float(cols[15]) if len(cols) > 15 else float("nan")
            except ValueError:
                bf_mean = float("nan")
            aa_len = max(0, (orf_len_nt - 3) // 3) if orf_len_nt > 0 else 0
            orf_class = classify_rpbp(orf_type)
            orf_class = reclassify_smorf(orf_class, aa_len)

            # Map Bayes factor (log scale, ~0 to 30+) to BED score 0-1000.
            if bf_mean == bf_mean:
                score_int = max(0, min(1000, int(round(min(bf_mean, 30.0) * 1000.0 / 30.0))))
            else:
                score_int = 0

            # Try to recover a transcript id from the BED 'name' field.
            tid = name.rsplit("_", 1)[0] if "_" in name else name
            gid = transcripts[tid].gene_id if tid in transcripts else ""

            orf_id = f"rpbp|{name}"
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
                    caller="rpbp",
                    sample_id=args.sample_id,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    gene_id=gid,
                    transcript_id=tid,
                    orf_class=orf_class,
                    aa_length=aa_len,
                    score=bf_mean if bf_mean == bf_mean else "",
                )
            )

    write_outputs(args.out_bed, args.out_tsv, bed_lines, tsv_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
