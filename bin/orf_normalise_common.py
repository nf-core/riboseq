#!/usr/bin/env python3
"""Shared utilities for ORF-normaliser scripts (ribotish/ribocode/ribotricer/rpbp).

Each per-caller script imports from this module. The module lives in `bin/`
which Nextflow adds to PATH; the scripts use `sys.path.insert(0, os.path.dirname(__file__))`
so the import works in the staged work directory.

Output BED12 column order:
    chrom start end name score strand thickStart thickEnd itemRgb
    blockCount blockSizes blockStarts

Output sidecar TSV columns:
    orf_id caller sample_id chrom start end strand gene_id transcript_id
    orf_class aa_length score
"""
from __future__ import annotations

import gzip
import io
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Iterator

BED12_HEADER = (
    "chrom\tstart\tend\tname\tscore\tstrand\t"
    "thickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts"
)
TSV_HEADER = (
    "orf_id\tcaller\tsample_id\tchrom\tstart\tend\tstrand\t"
    "gene_id\ttranscript_id\torf_class\taa_length\tscore"
)
# A simple attribute parser sufficient for GTF lines emitted by gffread,
# StringTie, Ensembl/GENCODE.
_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: list  # list of (start_0based, end_exclusive) genomic intervals, sorted ascending by start

    @property
    def length(self) -> int:
        return sum(e[1] - e[0] for e in self.exons)


def open_text(path: Path | str) -> IO[str]:
    """Open plain or gz text file."""
    p = Path(path)
    if str(p).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    return open(p, "r", encoding="utf-8")


def parse_attributes(field: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(field))


def load_transcripts(gtf_path: Path | str) -> dict[str, Transcript]:
    """Build a transcript-id -> Transcript map from a GTF, indexing exons only."""
    by_tid: dict[str, Transcript] = {}
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "exon":
                continue
            chrom = parts[0]
            start_1 = int(parts[3])
            end_1 = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            gid = attrs.get("gene_id", "")
            tx = by_tid.get(tid)
            if tx is None:
                tx = Transcript(
                    transcript_id=tid,
                    gene_id=gid,
                    chrom=chrom,
                    strand=strand,
                    exons=[],
                )
                by_tid[tid] = tx
            # convert 1-based inclusive GTF to 0-based half-open
            tx.exons.append((start_1 - 1, end_1))
    # Sort exons in genomic order (ascending start)
    for tx in by_tid.values():
        tx.exons.sort()
    return by_tid


def transcript_to_genomic_blocks(
    tx: Transcript, tx_start_0: int, tx_end_excl: int
) -> list[tuple[int, int]]:
    """Map a transcript-coordinate interval [tx_start_0, tx_end_excl)
    (in mRNA coordinates, 5'->3') onto genomic blocks.

    Returns a list of (genomic_start_0, genomic_end_excl) blocks, in genomic
    (ascending) order, ready for BED12 emission.
    """
    if tx_end_excl <= tx_start_0:
        return []
    length = tx.length
    if tx_start_0 < 0 or tx_end_excl > length:
        # Defensive: clamp to transcript bounds rather than fail loudly. ORF
        # callers occasionally emit a coordinate that crosses the polyA tail.
        tx_start_0 = max(0, tx_start_0)
        tx_end_excl = min(length, tx_end_excl)
        if tx_end_excl <= tx_start_0:
            return []

    # Build an exon walking order: 5'->3' on the mRNA. For '+' strand the
    # genomic order matches; for '-' strand we walk exons in reverse and the
    # mRNA-coord origin is the rightmost exon's end.
    blocks: list[tuple[int, int]] = []
    if tx.strand == "+":
        cumulative = 0
        for gs, ge in tx.exons:
            esize = ge - gs
            exon_mrna_start = cumulative
            exon_mrna_end = cumulative + esize
            # intersect mRNA-coord ORF with this exon
            lo = max(tx_start_0, exon_mrna_start)
            hi = min(tx_end_excl, exon_mrna_end)
            if hi > lo:
                gstart = gs + (lo - exon_mrna_start)
                gend = gs + (hi - exon_mrna_start)
                blocks.append((gstart, gend))
            cumulative = exon_mrna_end
    else:
        cumulative = 0
        # mRNA-coord origin is rightmost exon's end; walk exons high-to-low
        for gs, ge in reversed(tx.exons):
            esize = ge - gs
            exon_mrna_start = cumulative
            exon_mrna_end = cumulative + esize
            lo = max(tx_start_0, exon_mrna_start)
            hi = min(tx_end_excl, exon_mrna_end)
            if hi > lo:
                # On '-' strand, mRNA-coord offset `off` corresponds to genomic
                # position `ge - off - 1`. The mRNA-coord block [lo, hi) maps
                # to genomic [ge - (hi - exon_mrna_start), ge - (lo - exon_mrna_start))
                gend = ge - (lo - exon_mrna_start)
                gstart = ge - (hi - exon_mrna_start)
                blocks.append((gstart, gend))
            cumulative = exon_mrna_end
        blocks.sort()
    return blocks


def emit_bed12(
    chrom: str,
    blocks: list[tuple[int, int]],
    name: str,
    score: float | int | str,
    strand: str,
) -> str:
    """Format a BED12 line from genomic blocks."""
    if not blocks:
        return ""
    start = blocks[0][0]
    end = blocks[-1][1]
    block_count = len(blocks)
    block_sizes = ",".join(str(b[1] - b[0]) for b in blocks) + ","
    block_starts = ",".join(str(b[0] - start) for b in blocks) + ","
    # itemRgb '0' = unset
    return "\t".join(
        [
            chrom,
            str(start),
            str(end),
            name,
            str(score),
            strand,
            str(start),
            str(end),
            "0",
            str(block_count),
            block_sizes,
            block_starts,
        ]
    )


def emit_tsv_row(
    orf_id: str,
    caller: str,
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    gene_id: str,
    transcript_id: str,
    orf_class: str,
    aa_length: int | str,
    score: float | int | str,
) -> str:
    return "\t".join(
        [
            orf_id,
            caller,
            sample_id,
            chrom,
            str(start),
            str(end),
            strand,
            gene_id,
            transcript_id,
            orf_class,
            str(aa_length),
            str(score),
        ]
    )


def write_outputs(
    bed_path: Path | str,
    tsv_path: Path | str,
    bed_lines: Iterable[str],
    tsv_rows: Iterable[str],
) -> None:
    bed_lines = [l for l in bed_lines if l]
    with open(bed_path, "w") as bh:
        # No header on BED12 (downstream tools won't accept one).
        for l in bed_lines:
            bh.write(l + "\n")
    with open(tsv_path, "w") as th:
        th.write(TSV_HEADER + "\n")
        for r in tsv_rows:
            th.write(r + "\n")


def classify_ribotish(tis_type: str | None) -> str:
    """Map a Ribo-TISH `TisType` field to a unified orf_class label.

    Ribo-TISH TisType values seen in the wild include:
      'Annotated', 'CDSFrameOverlap', 'NCDSFrameOverlap', 'Internal',
      'Truncated:Annotated', 'Extended:Annotated', '5\'UTR', '3\'UTR',
      'Novel', 'Known' etc.
    """
    if not tis_type:
        return "other"
    t = tis_type.lower()
    if "5'utr" in t or "uorf" in t:
        return "uORF"
    if "3'utr" in t or "dorf" in t:
        return "dORF"
    if "annotated" in t and "truncated" not in t and "extended" not in t:
        return "canonical_cds"
    if "extended" in t or "truncated" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_ribocode(orf_type: str | None) -> str:
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "uorf" in t or "5'utr" in t:
        return "uORF"
    if "dorf" in t or "3'utr" in t:
        return "dORF"
    if "annotated" in t or "ccds" in t:
        return "canonical_cds"
    if "internal" in t:
        return "other"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_ribotricer(orf_type: str | None) -> str:
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "uorf" in t:
        return "uORF"
    if "dorf" in t:
        return "dORF"
    if "annotated" in t or "ccds" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_rpbp(orf_type: str | None) -> str:
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "five_prime" in t or "uorf" in t:
        return "uORF"
    if "three_prime" in t or "dorf" in t:
        return "dORF"
    if "canonical" in t or "annotated" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def reclassify_smorf(orf_class: str, aa_length: int) -> str:
    """Promote any ORF with aa_length <= 100 to smORF class, regardless of
    location-based classification."""
    if isinstance(aa_length, int) and aa_length <= 100 and aa_length > 0:
        return "smORF"
    return orf_class
