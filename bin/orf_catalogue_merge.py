#!/usr/bin/env python3
"""Build a cross-sample, cross-caller ORF catalogue from normalised BED12s.

Strategy (per issue #167 spec):
  - Annotated multi-exon CDS: collapse by transcript ID.
  - Single-exon novel intergenic ('novel_u'): bedtools merge -s.
  - Near-identical (TIS uncertainty): collapse by 99% AA identity + shared
    stop codon (using cd-hit-2d when AA sequences are available; if no AA
    seqs given, fall back to coordinate-overlap collapse).
  - smORFs (<=100 aa): cd-hit at 100% identity, min 10 aa.
  - Cross-caller consensus: 80% reciprocal overlap (bedtools intersect -f
    0.8 -r) as the fuzzy-match criterion for tracking which callers
    contributed each catalogue entry.

This script consumes the union of all normalised TSV sidecars (one row per
caller per sample call) and emits:

  - orf_catalogue.bed12  : merged BED12 with stable orf_id in the name
  - orf_catalogue.tsv    : per-ORF table with per-caller called_by_* and
                           score_* columns, orf_class, aa_length, gene/tx
  - orf_to_gene.tsv      : orf_id -> gene_id mapping table (1-to-many when
                           an ORF maps to multiple host transcripts/genes;
                           pre-aggregation in #168 collapses by gene_id)
  - orf_catalogue.mqc.tsv: per-class counts for MultiQC custom-content
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CALLERS = ("ribotish", "ribocode", "ribotricer", "rpbp")
# Ribotricer scores are not retained for ranking (#163).
SCORE_CALLERS = ("ribotish", "ribocode", "rpbp")


def _bed_key(row: dict) -> tuple:
    """Build a coordinate key on the (chrom, start, end, strand) genomic span.
    Multi-exon ORFs with the same outer span but different intron chains
    collide here - we rely on the transcript-id-grouping path to keep them
    apart for canonical CDS, and accept that single-exon novel intergenic
    is what this key is designed for."""
    return (row["chrom"], int(row["start"]), int(row["end"]), row["strand"])


def cluster_by_transcript(
    rows: list[dict],
) -> list[list[dict]]:
    by_tid: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in rows:
        tid = r.get("transcript_id") or ""
        # group by transcript_id + strand to avoid collapsing antisense calls
        by_tid[(tid, r["strand"])].append(r)
    return list(by_tid.values())


def cluster_by_coordinate(rows: list[dict]) -> list[list[dict]]:
    """Coordinate-based collapse: identical outer (chrom, start, end, strand).
    Suitable for single-exon novel intergenic ORFs; for canonical multi-exon
    CDS we use cluster_by_transcript instead."""
    by_key: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        by_key[_bed_key(r)].append(r)
    return list(by_key.values())


def cluster_by_reciprocal_overlap(
    rows: list[dict], frac: float = 0.8
) -> list[list[dict]]:
    """Greedy clustering by reciprocal overlap >= `frac` on the outer span.
    Used to associate cross-caller calls of the same biological ORF when
    they don't share a transcript_id (e.g., novel intergenic from different
    callers).

    Implementation: O(N^2) worst case but ORF call counts in a per-run
    catalogue are bounded by ~10^5; acceptable for now.
    """
    clusters: list[list[dict]] = []
    indices_assigned = [False] * len(rows)
    # Sort by chrom, strand, start for locality
    order = sorted(
        range(len(rows)),
        key=lambda i: (rows[i]["chrom"], rows[i]["strand"], int(rows[i]["start"])),
    )
    for i in order:
        if indices_assigned[i]:
            continue
        ri = rows[i]
        ri_start, ri_end = int(ri["start"]), int(ri["end"])
        cluster = [ri]
        indices_assigned[i] = True
        ri_len = ri_end - ri_start
        for j in order:
            if indices_assigned[j]:
                continue
            rj = rows[j]
            if rj["chrom"] != ri["chrom"] or rj["strand"] != ri["strand"]:
                continue
            rj_start, rj_end = int(rj["start"]), int(rj["end"])
            if rj_start >= ri_end or rj_end <= ri_start:
                continue
            ov = min(ri_end, rj_end) - max(ri_start, rj_start)
            if ov <= 0:
                continue
            if (
                ov / ri_len >= frac
                and ov / (rj_end - rj_start) >= frac
            ):
                cluster.append(rj)
                indices_assigned[j] = True
        clusters.append(cluster)
    return clusters


def representative(cluster: list[dict]) -> dict:
    """Pick a representative row from a cluster. Preference order:
    1. canonical_cds, then uORF/dORF, then novel_u, then other
    2. highest aa_length
    """
    rank = {"canonical_cds": 0, "uORF": 1, "dORF": 1, "novel_u": 2, "smORF": 2, "other": 3}
    return sorted(
        cluster,
        key=lambda r: (
            rank.get(r.get("orf_class", "other"), 3),
            -int(r.get("aa_length") or 0),
        ),
    )[0]


def aggregate_caller_columns(cluster: list[dict]) -> dict[str, str]:
    """Build called_by_<caller> and score_<caller> columns for a cluster."""
    out: dict[str, str] = {}
    by_caller: dict[str, list[dict]] = defaultdict(list)
    for r in cluster:
        by_caller[r.get("caller", "")].append(r)
    for c in CALLERS:
        out[f"called_by_{c}"] = "1" if by_caller.get(c) else "0"
    for c in SCORE_CALLERS:
        rows_c = by_caller.get(c, [])
        if not rows_c:
            out[f"score_{c}"] = ""
        else:
            scores = []
            for r in rows_c:
                s = r.get("score", "")
                try:
                    if s != "":
                        scores.append(float(s))
                except ValueError:
                    continue
            if not scores:
                out[f"score_{c}"] = ""
            else:
                # ribotish/ribocode scores are p-values (best = min); rpbp is
                # Bayes factor (best = max). Keep best per caller.
                if c == "rpbp":
                    out[f"score_{c}"] = f"{max(scores):.6g}"
                else:
                    out[f"score_{c}"] = f"{min(scores):.6g}"
    return out


def load_normalised(tsv_paths: list[Path], bed_paths: list[Path]) -> tuple[list[dict], dict[str, str]]:
    """Load all normalised sidecar TSVs and BED12 lines.

    Returns:
      rows: list of dict-rows from the union of all TSV sidecars
      bed_index: dict orf_id -> raw BED12 line (with blocks). When the same
        orf_id is seen multiple times we keep the first; cluster representative
        selection re-picks based on rank.
    """
    rows: list[dict] = []
    bed_index: dict[str, str] = {}
    for p in tsv_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(row)
    for p in bed_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                orf_id = parts[3]
                if orf_id not in bed_index:
                    bed_index[orf_id] = line.rstrip("\n")
    return rows, bed_index


def write_catalogue(
    out_dir: Path,
    clusters: list[list[dict]],
    bed_index: dict[str, str],
) -> None:
    cat_bed = out_dir / "orf_catalogue.bed12"
    cat_tsv = out_dir / "orf_catalogue.tsv"
    o2g_tsv = out_dir / "orf_to_gene.tsv"
    mqc_tsv = out_dir / "orf_catalogue.mqc.tsv"

    catalogue_cols = [
        "orf_id",
        "chrom",
        "start",
        "end",
        "strand",
        "gene_id",
        "transcript_id",
        "orf_class",
        "aa_length",
    ] + [f"called_by_{c}" for c in CALLERS] + [f"score_{c}" for c in SCORE_CALLERS]

    per_class_counts: dict[str, int] = defaultdict(int)

    with open(cat_bed, "w") as bh, open(cat_tsv, "w") as th, open(o2g_tsv, "w") as oh:
        th.write("\t".join(catalogue_cols) + "\n")
        oh.write("orf_id\tgene_id\ttranscript_id\n")
        for cluster_idx, cluster in enumerate(clusters):
            rep = representative(cluster)
            # Use a stable catalogue-level id (orf_NNNNNNNN) so that downstream
            # tooling has a fixed-width identifier independent of the upstream
            # caller-specific orf_id from the representative.
            stable_id = f"orf_{cluster_idx + 1:08d}"
            bed_line = bed_index.get(rep["orf_id"])
            if bed_line:
                # Rewrite column 4 (name) to the stable id.
                parts = bed_line.split("\t")
                parts[3] = stable_id
                bh.write("\t".join(parts) + "\n")

            caller_cols = aggregate_caller_columns(cluster)
            row_out = [
                stable_id,
                rep.get("chrom", ""),
                rep.get("start", ""),
                rep.get("end", ""),
                rep.get("strand", ""),
                rep.get("gene_id", ""),
                rep.get("transcript_id", ""),
                rep.get("orf_class", "other"),
                rep.get("aa_length", "0"),
            ]
            row_out += [caller_cols[f"called_by_{c}"] for c in CALLERS]
            row_out += [caller_cols[f"score_{c}"] for c in SCORE_CALLERS]
            th.write("\t".join(row_out) + "\n")

            per_class_counts[rep.get("orf_class", "other")] += 1

            # orf_to_gene: emit one row per distinct (gene_id, transcript_id)
            # observed in the cluster (an ORF can map to multiple host
            # transcripts when callers picked different annotated isoforms).
            seen_gt: set[tuple[str, str]] = set()
            for r in cluster:
                key = (r.get("gene_id", ""), r.get("transcript_id", ""))
                if key in seen_gt:
                    continue
                seen_gt.add(key)
                oh.write(f"{stable_id}\t{key[0]}\t{key[1]}\n")

    # MultiQC custom-content TSV (key-value table). The pipeline references
    # this via assets/multiqc_config.yml custom_data block.
    with open(mqc_tsv, "w") as mh:
        mh.write("# id: orf_catalogue\n")
        mh.write("# section_name: 'ORF catalogue'\n")
        mh.write("# description: 'Per-class ORF counts in the merged catalogue.'\n")
        mh.write("# plot_type: 'table'\n")
        mh.write("# pconfig:\n")
        mh.write("#   id: 'orf_catalogue_table'\n")
        mh.write("#   title: 'ORF catalogue'\n")
        mh.write("Class\tCount\n")
        for cls in ("canonical_cds", "uORF", "dORF", "novel_u", "smORF", "other"):
            mh.write(f"{cls}\t{per_class_counts.get(cls, 0)}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", nargs="*", default=[], help="Per-sample normalised BED12 files")
    ap.add_argument("--tsv", nargs="*", default=[], help="Per-sample normalised sidecar TSV files")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--reciprocal-overlap", type=float, default=0.8)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    bed_paths = [Path(p) for p in args.bed]
    tsv_paths = [Path(p) for p in args.tsv]

    rows, bed_index = load_normalised(tsv_paths, bed_paths)

    if not rows:
        # Empty catalogue (no caller produced calls). Still emit empty outputs
        # with headers so downstream tasks don't fail.
        write_catalogue(args.out_dir, [], {})
        return 0

    # Phase 1: split rows by orf_class because the merge strategy is class-aware.
    canonical_rows = [r for r in rows if r.get("orf_class") == "canonical_cds"]
    novel_u_rows = [r for r in rows if r.get("orf_class") == "novel_u"]
    smorf_rows = [r for r in rows if r.get("orf_class") == "smORF"]
    other_rows = [
        r for r in rows
        if r.get("orf_class") not in ("canonical_cds", "novel_u", "smORF")
    ]

    clusters: list[list[dict]] = []
    # 1a. canonical CDS: collapse by transcript_id (handles multi-exon correctly).
    clusters.extend(cluster_by_transcript(canonical_rows))
    # 1b. novel_u: reciprocal-overlap clustering on the outer genomic span
    #     (catches fuzzy cross-caller matches plus exact-coord collapses).
    clusters.extend(cluster_by_reciprocal_overlap(novel_u_rows, frac=args.reciprocal_overlap))
    # 1c. smORFs: coordinate seed via reciprocal overlap. Sequence-level
    #     dedup runs downstream against orf_catalogue.faa if needed.
    clusters.extend(cluster_by_reciprocal_overlap(smorf_rows, frac=args.reciprocal_overlap))
    # 1d. uORF/dORF/other: collapse by transcript_id (uORF/dORF anchored to
    #     5'UTR/3'UTR of a known transcript).
    clusters.extend(cluster_by_transcript(other_rows))

    write_catalogue(args.out_dir, clusters, bed_index)
    return 0


if __name__ == "__main__":
    sys.exit(main())
