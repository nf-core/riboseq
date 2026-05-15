#!/usr/bin/env python3
"""Normalise a Ribo-TISH predict output to unified BED12 + sidecar TSV.

Ribo-TISH `predict` writes a tab-separated table with a header row. Relevant
columns (header names; case as emitted by ribotish 0.2.x):

    Tid          - transcript id
    Gid          - gene id (only when -a/--gtf2 is set)
    GenomePos    - 'chr:start-end:strand' (1-based inclusive, genome coords)
    StartCodon   - e.g. 'ATG', 'CTG'
    Codon        - e.g. 'ATG', 'CTG'
    AALen        - amino acid length
    TisType      - classification (Annotated, 5'UTR, 3'UTR, Novel, Internal, ...)
    Blocks       - 'a:b,c:d,...' (1-based inclusive) when ribotish emits blocks
    InFrameCount - in-frame P-site count
    FPKM         - FPKM
    Pvalue       - global combined p-value
    Pvalcombined - combined p-value (newer ribotish)

The script prefers the `Blocks` column for BED12 construction; when not
present (older ribotish output), it falls back to looking up the transcript
in the GTF and intersecting the ORF transcript-coord interval with its exons.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_normalise_common import (
    classify_ribotish,
    emit_bed12,
    emit_tsv_row,
    load_transcripts,
    reclassify_smorf,
    transcript_to_genomic_blocks,
    write_outputs,
)

_GENPOS_RE = re.compile(r"^([^:]+):(\d+)-(\d+):([+\-])$")


def parse_genomepos(s: str) -> tuple[str, int, int, str] | None:
    if not s:
        return None
    m = _GENPOS_RE.match(s.strip())
    if not m:
        return None
    chrom, start_1, end_1, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    # 1-based inclusive -> 0-based half-open
    return chrom, start_1 - 1, end_1, strand


def parse_blocks(s: str) -> list[tuple[int, int]]:
    """Parse Ribo-TISH 'Blocks' field. Format: 'a:b,c:d,...' 1-based inclusive.

    Returns list of (start_0, end_excl) genomic intervals, sorted ascending.
    """
    if not s or s in ("-", "."):
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            a, b = tok.split(":", 1)
        elif "-" in tok:
            a, b = tok.split("-", 1)
        else:
            continue
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

    # If the input is empty (zero-call sample, stub run), emit empty outputs.
    if not args.input.exists() or args.input.stat().st_size == 0:
        write_outputs(args.out_bed, args.out_tsv, [], [])
        return 0

    transcripts = load_transcripts(args.gtf)

    bed_lines: list[str] = []
    tsv_rows: list[str] = []
    seen_orf_ids: set[str] = set()

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid = row.get("Tid", "") or ""
            gid = row.get("Gid", "") or (transcripts.get(tid).gene_id if tid in transcripts else "")
            try:
                aa_len = int(row.get("AALen", "0") or 0)
            except ValueError:
                aa_len = 0
            tis_type = row.get("TisType", "") or row.get("TISType", "")
            orf_class = classify_ribotish(tis_type)
            orf_class = reclassify_smorf(orf_class, aa_len)

            # Prefer combined p-value as a score proxy; smaller = better.
            score_raw = row.get("Pvalcombined") or row.get("Pvalue") or row.get("Pvalcom") or "1"
            try:
                pval = float(score_raw)
                # BED 'score' int 0-1000; map 1-Pval to 0-1000.
                score_int = max(0, min(1000, int(round((1.0 - pval) * 1000))))
            except (TypeError, ValueError):
                pval = float("nan")
                score_int = 0

            gp = parse_genomepos(row.get("GenomePos", ""))
            if gp is None:
                continue
            chrom, start, end, strand = gp

            blocks = parse_blocks(row.get("Blocks", ""))
            if not blocks:
                # Fall back to exon-walking against the GTF if Blocks absent.
                tx = transcripts.get(tid)
                if tx is None:
                    blocks = [(start, end)]
                else:
                    # Map the ORF's transcript-coord range to genomic blocks. We
                    # don't know the transcript-coord offset directly; reconstruct
                    # by intersecting the ORF genomic [start, end) with the tx
                    # exons - this gives correct multi-exon blocks for canonical
                    # CDS ORFs, and a single block for novel intergenic.
                    blocks = []
                    for gs, ge in tx.exons:
                        lo = max(start, gs)
                        hi = min(end, ge)
                        if hi > lo:
                            blocks.append((lo, hi))
                    if not blocks:
                        blocks = [(start, end)]

            orf_id = f"ribotish|{tid}|{chrom}:{blocks[0][0]}-{blocks[-1][1]}:{strand}"
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
                    caller="ribotish",
                    sample_id=args.sample_id,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    gene_id=gid,
                    transcript_id=tid,
                    orf_class=orf_class,
                    aa_length=aa_len,
                    score=pval if pval == pval else "",  # NaN -> empty
                )
            )

    write_outputs(args.out_bed, args.out_tsv, bed_lines, tsv_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
