#!/usr/bin/env python3
"""Normalise a RiboCode per-sample ORF table to unified BED12 + sidecar TSV.

RiboCode writes a tab-separated table. Relevant columns (RiboCode 1.2.x):

    ORF_ID            - composite id (transcript_id_ORF_start_ORF_stop)
    ORF_type          - annotated|uORF|dORF|Overlap_uORF|Overlap_dORF|Internal|...
    transcript_id     - source transcript id
    gene_id           - source gene id
    chrom             - chromosome
    strand            - +/-
    ORF_tstart        - 1-based transcript start (5' end of ORF)
    ORF_tstop         - 1-based transcript stop (3' end, inclusive)
    ORF_length        - nucleotide length (3 * aa + 3 for stop)
    AA_length         - amino acid length (excluding stop)
    pval_combined     - combined p-value across read-length bins (lower = better)
    pval_frame0/1/2   - per-frame p-values

The script uses ORF_tstart/ORF_tstop with the transcript's exon model to
produce genomic BED12 blocks (handles multi-exon canonical CDS correctly).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_normalise_common import (
    classify_ribocode,
    emit_bed12,
    emit_tsv_row,
    load_transcripts,
    reclassify_smorf,
    transcript_to_genomic_blocks,
    write_outputs,
)


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
            tid = row.get("transcript_id") or row.get("Transcript_id") or ""
            gid = row.get("gene_id") or row.get("Gene_id") or ""
            chrom = row.get("chrom") or row.get("Chrom") or ""
            strand = row.get("strand") or row.get("Strand") or "+"
            try:
                t_start_1 = int(row.get("ORF_tstart") or row.get("ORF_Tstart") or 0)
                t_stop_1 = int(row.get("ORF_tstop") or row.get("ORF_Tstop") or 0)
            except ValueError:
                continue
            try:
                aa_len = int(row.get("AA_length") or row.get("AAlength") or 0)
            except ValueError:
                aa_len = 0
            orf_type = row.get("ORF_type") or row.get("Type") or ""

            orf_class = classify_ribocode(orf_type)
            orf_class = reclassify_smorf(orf_class, aa_len)

            score_raw = row.get("pval_combined") or row.get("Pval_combined") or "1"
            try:
                pval = float(score_raw)
                score_int = max(0, min(1000, int(round((1.0 - pval) * 1000))))
            except (TypeError, ValueError):
                pval = float("nan")
                score_int = 0

            tx = transcripts.get(tid)
            blocks: list[tuple[int, int]] = []
            if tx is not None and t_stop_1 >= t_start_1 > 0:
                # ORF_tstart/ORF_tstop are 1-based inclusive transcript coords
                blocks = transcript_to_genomic_blocks(
                    tx, t_start_1 - 1, t_stop_1
                )
                if not chrom:
                    chrom = tx.chrom
                if not strand or strand not in ("+", "-"):
                    strand = tx.strand
            if not blocks:
                # Fall back to the ORF_ID's genomic coordinates if RiboCode
                # supplied genome columns directly.
                g_start = row.get("ORF_gstart") or row.get("ORF_genomic_start") or ""
                g_stop = row.get("ORF_gstop") or row.get("ORF_genomic_stop") or ""
                try:
                    g_start_i = int(g_start)
                    g_stop_i = int(g_stop)
                    if g_stop_i < g_start_i:
                        g_start_i, g_stop_i = g_stop_i, g_start_i
                    blocks = [(g_start_i - 1, g_stop_i)]
                except (TypeError, ValueError):
                    continue

            orf_id = row.get("ORF_ID") or f"ribocode|{tid}|{t_start_1}-{t_stop_1}"
            orf_id = f"ribocode|{orf_id}"
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
                    caller="ribocode",
                    sample_id=args.sample_id,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    gene_id=gid,
                    transcript_id=tid,
                    orf_class=orf_class,
                    aa_length=aa_len,
                    score=pval if pval == pval else "",
                )
            )

    write_outputs(args.out_bed, args.out_tsv, bed_lines, tsv_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
