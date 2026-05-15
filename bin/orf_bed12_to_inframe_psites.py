#!/usr/bin/env python3
"""Convert a BED12 ORF catalogue to a BED6 of in-frame P-site positions.

For each ORF (BED12 record), this script emits one BED6 row per codon-start
position along the mRNA-spliced ORF length, projected back onto genomic
coordinates. Frame is defined by the ORF's own start codon (ATG = frame 0),
independent of any GTF `phase` annotation, because catalogue ORFs may
originate from novel transcripts or non-canonical starts where GTF phase is
unset or unreliable.

BED12 conventions assumed (UCSC):
  chrom start end name score strand thickStart thickEnd itemRgb
  blockCount blockSizes blockStarts
where `start` is 0-based and `end` is half-open. blockStarts are offsets
from `start`. Blocks are listed in ascending genomic-coordinate order on
both strands; mRNA-order traversal is therefore left-to-right for `+`
strand ORFs and right-to-left for `-` strand ORFs.

Output (BED6, 0-based half-open):
  chrom  start  end  orf_id  0  strand

Each row is a single-nucleotide interval at the A position of an in-frame
codon (i.e. the 5' end of a codon on the mRNA). Duplicate rows are
suppressed (an ORF cannot claim the same in-frame position twice).
"""
from __future__ import annotations

import argparse
import sys


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bed12", help="Input BED12 ORF catalogue (use '-' for stdin).")
    p.add_argument("-o", "--output", default="-", help="Output BED6 path (use '-' for stdout).")
    return p.parse_args(argv)


def expand_blocks(start: int, block_sizes: list[int], block_starts: list[int]) -> list[tuple[int, int]]:
    """Return list of (block_start, block_end) genomic intervals, half-open, ascending."""
    blocks = []
    for sz, off in zip(block_sizes, block_starts):
        bs = start + off
        be = bs + sz
        blocks.append((bs, be))
    blocks.sort()
    return blocks


def codon_start_positions(blocks: list[tuple[int, int]], strand: str) -> list[int]:
    """Yield genomic 0-based positions for each in-frame codon-start nucleotide.

    Walks the spliced ORF in mRNA order (5' -> 3'), emitting every 3rd
    nucleotide starting from mRNA position 0. For each emitted mRNA
    position, computes the corresponding genomic coordinate.
    """
    # Pre-compute (block_index, genomic_start, genomic_end, mrna_offset_at_block_start)
    # in mRNA traversal order.
    if strand == "+":
        ordered_blocks = list(blocks)
    elif strand == "-":
        ordered_blocks = list(reversed(blocks))
    else:
        return []

    total_len = sum(b[1] - b[0] for b in ordered_blocks)
    positions: list[int] = []

    # Walk by codon (step of 3 along mRNA). For each codon-start mRNA index,
    # find which block it lies in and convert to genomic coordinate.
    cumulative = 0
    block_iter_idx = 0
    blk_start, blk_end = ordered_blocks[0]
    blk_len = blk_end - blk_start
    blk_cum_start = 0

    for mrna_pos in range(0, total_len, 3):
        # Advance to the block containing mrna_pos.
        while mrna_pos >= blk_cum_start + blk_len:
            block_iter_idx += 1
            if block_iter_idx >= len(ordered_blocks):
                return positions
            blk_cum_start += blk_len
            blk_start, blk_end = ordered_blocks[block_iter_idx]
            blk_len = blk_end - blk_start

        offset_in_block = mrna_pos - blk_cum_start
        if strand == "+":
            genomic = blk_start + offset_in_block
        else:
            # mRNA 5' end of the block on '-' strand maps to the block's
            # genomic end - 1.
            genomic = blk_end - 1 - offset_in_block
        positions.append(genomic)

    return positions


def process(line: str) -> list[str]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 12:
        return []
    chrom = parts[0]
    start = int(parts[1])
    # end = int(parts[2])  # unused; reconstructed from blocks
    name = parts[3]
    strand = parts[5]
    block_count = int(parts[9])
    block_sizes = [int(x) for x in parts[10].rstrip(",").split(",") if x != ""]
    block_starts = [int(x) for x in parts[11].rstrip(",").split(",") if x != ""]
    if len(block_sizes) != block_count or len(block_starts) != block_count:
        # Malformed record; skip silently. The catalogue producer should
        # have caught this.
        return []
    blocks = expand_blocks(start, block_sizes, block_starts)
    out_lines: list[str] = []
    seen: set[int] = set()
    for gpos in codon_start_positions(blocks, strand):
        if gpos in seen:
            continue
        seen.add(gpos)
        out_lines.append(f"{chrom}\t{gpos}\t{gpos + 1}\t{name}\t0\t{strand}")
    return out_lines


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    in_fh = sys.stdin if args.bed12 == "-" else open(args.bed12)
    out_fh = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        for line in in_fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            for row in process(line):
                out_fh.write(row + "\n")
    finally:
        if in_fh is not sys.stdin:
            in_fh.close()
        if out_fh is not sys.stdout:
            out_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
