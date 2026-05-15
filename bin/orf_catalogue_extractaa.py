#!/usr/bin/env python3
"""Extract AA sequences for all ORFs in a unified BED12 catalogue.

Inputs:
  --bed12   : merged catalogue (BED12 with stable orf_NNN ids in column 4)
  --fasta   : genome FASTA (indexed; we'll create the .fai if missing)

Output:
  --out-faa : FASTA of one record per ORF (>orf_NNN <aa_length=N>; sequence)

Translation uses the standard genetic code. ORF blocks are concatenated in
mRNA-coord order (5'->3') before translation. The starting block depends on
strand: '+' uses BED12 blocks left-to-right, '-' reverses and reverse-
complements before translation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(nt: str) -> str:
    nt = nt.upper()
    aa = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i:i+3]
        if len(codon) < 3:
            break
        aa.append(CODON_TABLE.get(codon, "X"))
    # Drop trailing stop codon for AA length consistency with caller outputs.
    if aa and aa[-1] == "*":
        aa = aa[:-1]
    return "".join(aa)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed12", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--out-faa", required=True, type=Path)
    args = ap.parse_args()

    if not args.bed12.exists() or args.bed12.stat().st_size == 0:
        args.out_faa.write_text("")
        return 0

    fa = pysam.FastaFile(str(args.fasta))

    with open(args.out_faa, "w") as oh:
        with open(args.bed12) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                chrom = parts[0]
                start = int(parts[1])
                name = parts[3]
                strand = parts[5]
                try:
                    block_count = int(parts[9])
                except ValueError:
                    continue
                block_sizes = [int(x) for x in parts[10].rstrip(",").split(",") if x]
                block_starts = [int(x) for x in parts[11].rstrip(",").split(",") if x]
                if len(block_sizes) != block_count or len(block_starts) != block_count:
                    continue

                nt_parts: list[str] = []
                for bs, sz in zip(block_starts, block_sizes):
                    s = start + bs
                    e = s + sz
                    try:
                        seq = fa.fetch(chrom, s, e)
                    except (KeyError, ValueError):
                        seq = ""
                    nt_parts.append(seq)
                if strand == "-":
                    nt = revcomp("".join(nt_parts))
                else:
                    nt = "".join(nt_parts)
                aa = translate(nt)
                if not aa:
                    continue
                oh.write(f">{name} aa_length={len(aa)}\n{aa}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
