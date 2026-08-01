#!/usr/bin/env python3
"""Merge per-sample, per-caller normalised ORF BED12s into a unified catalogue.

Rows are partitioned by clustering *strategy*. The harmonised `orf_class`
from custom/orfnormalise selects the strategy but is never part of a
grouping key, because callers disagree on class for the same ORF and keying
on it would emit one catalogue row per disagreeing caller:

  - canonical_cds: grouped by (transcript_id, strand), then reciprocal-overlap
    clustered within the transcript. Multi-exon CDSs from different callers
    still meet, while a short truncated variant stays a separate row instead
    of being folded into the full-length CDS.
  - uORF, uoORF, dORF, doORF, intORF, other: collapse by (transcript_id,
    strand, start, end). A transcript can host several, so the outer span
    joins the key.
  - novel_u: greedy reciprocal-overlap clustering on the outer genomic span
    at `--reciprocal-overlap` (default 0.8). Catches fuzzy cross-caller
    matches and exact-coordinate collapses together.

Every input row must land in exactly one cluster; the run aborts if the
class vocabulary grows beyond CLASS_ORDER or if any row goes unassigned.

Cross-caller consensus is recorded in two column families on the output
catalogue TSV:

  - `called_by_<caller>`: 0/1 indicator per supported caller.
  - `score_<caller>`:     best score from that caller within the cluster.
                          Score direction is per-caller (p-values are
                          minimised; Bayes factors / phase scores are
                          maximised).

Cross-sample recurrence is recorded in two further columns:

  - `n_samples`: number of distinct samples contributing to the cluster.
  - `samples`:   sorted, comma-separated list of those sample ids.

Outputs:

  ${prefix}.bed12                merged catalogue (genomic blocks).
  ${prefix}.tsv                  per-ORF table with caller-tracking cols.
  ${prefix}.orf_to_gene.tsv      one row per (orf_id, gene_id, transcript_id);
                                 an ORF can map to multiple host transcripts.
  ${prefix}.mqc.tsv              MultiQC custom-content sidecar
                                 (per-class counts).
"""

import argparse
import csv
import glob
import platform
import shlex
import sys
from collections import defaultdict
from pathlib import Path

import yaml

CALLERS = ("ribotish", "ribocode", "ribotricer", "rpbp", "price")

# Direction of each caller's native score: 'min' = lower is better
# (p-values); 'max' = higher is better (Bayes factors, phase scores).
# When aggregating per-cluster scores in the catalogue TSV, we keep the
# direction-appropriate best.
SCORE_DIRECTIONS = {
    "ribotish": "min",
    "ribocode": "min",
    "ribotricer": "max",
    "rpbp": "max",
    "price": "min",
}

# Display order for the MultiQC per-class table.
CLASS_ORDER = ("canonical_cds", "uORF", "uoORF", "dORF", "doORF", "intORF", "novel_u", "other")

# Representative preference when a cluster's members disagree on class, most
# specific first: a caller that resolved the CDS-overlap wins over one that
# could not. Must cover exactly CLASS_ORDER; asserted at runtime.
CLASS_SPECIFICITY = ("canonical_cds", "uoORF", "uORF", "doORF", "dORF", "intORF", "novel_u", "other")

# Clustering strategy per class. `orf_class` is deliberately not part of any
# grouping key: callers disagree on class for the same ORF (Ribo-TISH cannot
# report the CDS-overlapping uORF form at all), so keying on it would split
# those calls into separate catalogue rows.
OVERLAP_CLASSES = ("novel_u",)
LOOSE_CLASSES = ("canonical_cds",)


def group_by(rows, keyfn):
    """Group rows into clusters keyed by ``keyfn(row)``, preserving first-seen order."""
    groups = defaultdict(list)
    for r in rows:
        groups[keyfn(r)].append(r)
    return list(groups.values())


def cluster_by_reciprocal_overlap(rows, frac=0.8):
    """Greedy clustering by reciprocal overlap >= ``frac`` on the outer span.

    Used to merge cross-caller calls of the same biological ORF when they
    don't share a transcript_id (e.g. novel intergenic from different
    callers). O(N^2) worst case but bounded by per-run ORF counts.

    Greedy and order-dependent: at the default frac=0.8 a chain A-B-C
    where A overlaps B at 0.85, B overlaps C at 0.85, but A overlaps C
    at 0.75 can either land as {A, B, C} or {A, B} + {C} depending on
    iteration order. Acceptable in practice (overlap chains that
    straddle the threshold are rare at frac=0.8) but worth flagging for
    consumers that want strictly transitive clustering.
    """
    clusters = []
    assigned = [False] * len(rows)
    order = sorted(
        range(len(rows)),
        key=lambda i: (rows[i]["chrom"], rows[i]["strand"], int(rows[i]["start"])),
    )
    for i in order:
        if assigned[i]:
            continue
        ri = rows[i]
        ri_start, ri_end = int(ri["start"]), int(ri["end"])
        ri_len = ri_end - ri_start
        cluster = [ri]
        assigned[i] = True
        for j in order:
            if assigned[j]:
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
            rj_len = rj_end - rj_start
            if ri_len > 0 and rj_len > 0 and ov / ri_len >= frac and ov / rj_len >= frac:
                cluster.append(rj)
                assigned[j] = True
        clusters.append(cluster)
    return clusters


def representative(cluster):
    """Pick a representative row from a cluster.

    Class preference follows CLASS_SPECIFICITY; ties broken by longest aa_length.
    """
    rank = {c: i for i, c in enumerate(CLASS_SPECIFICITY)}
    return sorted(
        cluster,
        key=lambda r: (rank.get(r.get("orf_class", "other"), len(rank)), -int(r.get("aa_length") or 0)),
    )[0]


def aggregate_caller_columns(cluster):
    out = {}
    by_caller = defaultdict(list)
    for r in cluster:
        by_caller[r.get("caller", "")].append(r)
    for c in CALLERS:
        out[f"called_by_{c}"] = "1" if by_caller.get(c) else "0"
    for c in CALLERS:
        rows_c = by_caller.get(c, [])
        if not rows_c:
            out[f"score_{c}"] = ""
            continue
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
            best = max(scores) if SCORE_DIRECTIONS[c] == "max" else min(scores)
            out[f"score_{c}"] = f"{best:.6g}"
    return out


def load_normalised(tsv_paths, bed_paths):
    rows = []
    bed_index = {}
    for p in tsv_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p, newline="") as fh:
            data = (line for line in fh if not line.startswith("#"))
            rows.extend(csv.DictReader(data, delimiter="\\t"))
    for p in bed_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p) as fh:
            for line in fh:
                parts = line.rstrip("\\n").split("\\t")
                if len(parts) < 12:
                    continue
                orf_id = parts[3]
                if orf_id not in bed_index:
                    bed_index[orf_id] = line.rstrip("\\n")
    return rows, bed_index


def write_catalogue(prefix, clusters, bed_index, min_callers=1, min_samples=1):
    cat_bed = Path(f"{prefix}.bed12")
    cat_tsv = Path(f"{prefix}.tsv")
    o2g_tsv = Path(f"{prefix}.orf_to_gene.tsv")
    mqc_tsv = Path(f"{prefix}.mqc.tsv")
    # Consensus view: the same catalogue filtered to ORFs supported by at least
    # `min_callers` distinct callers and recurring in at least `min_samples`
    # samples. At the defaults (1, 1) it is identical to the full catalogue.
    con_bed = Path(f"{prefix}.consensus.bed12")
    con_tsv = Path(f"{prefix}.consensus.tsv")
    con_o2g = Path(f"{prefix}.consensus.orf_to_gene.tsv")

    catalogue_cols = (
        ["orf_id", "chrom", "start", "end", "strand", "gene_id", "transcript_id", "orf_class", "aa_length"]
        + [f"called_by_{c}" for c in CALLERS]
        + [f"score_{c}" for c in CALLERS]
        + ["n_samples", "samples", "orf_type_native", "is_smorf"]
    )

    per_class_counts = defaultdict(int)

    # orf_ids are assigned in iteration order below, so sort first to make the
    # numbering deterministic.
    def _sort_key(cluster):
        r = representative(cluster)
        return (
            r.get("chrom", ""),
            int(r.get("start") or 0),
            int(r.get("end") or 0),
            r.get("strand", ""),
            r.get("gene_id") or "",
            r.get("transcript_id") or "",
            r.get("orf_class", ""),
        )

    clusters = sorted(clusters, key=_sort_key)

    header = "\\t".join(catalogue_cols) + "\\n"
    o2g_header = "orf_id\\tgene_id\\ttranscript_id\\n"

    # First pass: write the full catalogue, recording which rows qualify for the
    # consensus view (>= min_callers distinct callers and >= min_samples samples).
    consensus_rows = []
    with open(cat_bed, "w") as bh, open(cat_tsv, "w") as th, open(o2g_tsv, "w") as oh:
        th.write(header)
        oh.write(o2g_header)
        for idx, cluster in enumerate(clusters):
            rep = representative(cluster)
            stable_id = f"orf_{idx + 1:08d}"
            bed_out = None
            bed_line = bed_index.get(rep["orf_id"])
            if bed_line:
                parts = bed_line.split("\\t")
                parts[3] = stable_id
                bed_out = "\\t".join(parts) + "\\n"
                bh.write(bed_out)

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
            row_out += [caller_cols[f"score_{c}"] for c in CALLERS]

            sample_ids = sorted({r.get("sample_id", "") for r in cluster if r.get("sample_id")})
            row_out += [str(len(sample_ids)), ",".join(sample_ids)]
            # Every native label in the cluster, so a cross-caller class
            # disagreement stays auditable after the harmonised class is picked.
            natives = sorted({(r.get("orf_type_native") or "").strip() for r in cluster} - {""})
            row_out += [",".join(natives)]
            # Carried from the representative, so is_smorf always agrees with
            # the aa_length emitted on the same row.
            row_out += [rep.get("is_smorf", "0")]
            row_line = "\\t".join(row_out) + "\\n"
            th.write(row_line)

            per_class_counts[rep.get("orf_class", "other")] += 1

            o2g_pairs = []
            seen_gt = set()
            for r in cluster:
                key = (r.get("gene_id", ""), r.get("transcript_id", ""))
                if key in seen_gt:
                    continue
                seen_gt.add(key)
                oh.write(f"{stable_id}\\t{key[0]}\\t{key[1]}\\n")
                o2g_pairs.append(key)

            # Consensus membership: distinct callers that flagged the cluster
            # and distinct contributing samples must each meet the threshold.
            n_callers = sum(1 for c in CALLERS if caller_cols[f"called_by_{c}"] == "1")
            if n_callers >= min_callers and len(sample_ids) >= min_samples:
                consensus_rows.append((bed_out, row_line, stable_id, o2g_pairs))

    # Second pass: the consensus view is the full catalogue restricted to the
    # qualifying rows.
    with open(con_bed, "w") as cbh, open(con_tsv, "w") as cth, open(con_o2g, "w") as coh:
        cth.write(header)
        coh.write(o2g_header)
        for bed_out, row_line, stable_id, o2g_pairs in consensus_rows:
            if bed_out:
                cbh.write(bed_out)
            cth.write(row_line)
            for gid, tid in o2g_pairs:
                coh.write(f"{stable_id}\\t{gid}\\t{tid}\\n")

    with open(mqc_tsv, "w") as mh:
        mh.write("# id: orf_catalogue\\n")
        mh.write("# section_name: 'ORF catalogue'\\n")
        mh.write("# description: 'Per-class ORF counts in the merged catalogue.'\\n")
        mh.write("# plot_type: 'table'\\n")
        mh.write("# pconfig:\\n")
        mh.write("#   id: 'orf_catalogue_table'\\n")
        mh.write("#   title: 'ORF catalogue'\\n")
        mh.write("Class\\tCount\\n")
        for cls in CLASS_ORDER:
            mh.write(f"{cls}\\t{per_class_counts.get(cls, 0)}\\n")


def write_versions():
    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {"${task.process}": {"python": platform.python_version()}},
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        default=0.8,
        help="Reciprocal-overlap fraction for novel_u and within-transcript CDS clustering (default: 0.8)",
    )
    parser.add_argument(
        "--min-callers",
        type=int,
        default=1,
        help="Minimum distinct callers for an ORF to enter the consensus view (default: 1)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum distinct samples for an ORF to enter the consensus view (default: 1)",
    )
    args = parser.parse_args(shlex.split("${args}"))

    if set(CLASS_SPECIFICITY) != set(CLASS_ORDER):
        sys.exit("orfmerge: CLASS_SPECIFICITY and CLASS_ORDER must cover the same classes")

    bed_paths = sorted(Path(p) for p in glob.glob("beds/*"))
    tsv_paths = sorted(Path(p) for p in glob.glob("tsvs/*"))

    prefix = "${prefix}"
    rows, bed_index = load_normalised(tsv_paths, bed_paths)

    if not rows:
        write_catalogue(prefix, [], {}, min_callers=args.min_callers, min_samples=args.min_samples)
        write_versions()
        return 0

    by_class = defaultdict(list)
    for r in rows:
        by_class[r.get("orf_class", "other")].append(r)

    unknown = sorted(set(by_class) - set(CLASS_ORDER))
    if unknown:
        sys.exit(f"orfmerge: unknown orf_class value(s) {unknown}; update CLASS_ORDER")

    clusters = []
    # Not transcript-anchored: one reciprocal-overlap pass over the union of
    # these classes, so calls disagreeing on class still reach one cluster.
    non_anchored = [r for cls in OVERLAP_CLASSES for r in by_class.get(cls, [])]
    clusters.extend(cluster_by_reciprocal_overlap(non_anchored, frac=args.reciprocal_overlap))
    # Grouped by transcript, then overlap-clustered within the transcript.
    # Grouping on (transcript_id, strand) alone would fold a short truncated
    # CDS variant into the full-length CDS and emit only the longest, so the
    # overlap pass keeps genuinely different ORFs on one transcript apart
    # while still merging cross-caller calls whose bounds differ slightly.
    loose = [r for cls in LOOSE_CLASSES for r in by_class.get(cls, [])]
    for grp in group_by(loose, lambda r: (r.get("transcript_id") or "", r["strand"])):
        clusters.extend(cluster_by_reciprocal_overlap(grp, frac=args.reciprocal_overlap))
    # Transcript-anchored, but a transcript can host several, so the outer
    # span joins the key.
    keyed = set(OVERLAP_CLASSES) | set(LOOSE_CLASSES)
    anchored = [r for cls, rows_c in by_class.items() if cls not in keyed for r in rows_c]
    clusters.extend(
        group_by(
            anchored,
            lambda r: (r.get("transcript_id") or "", r["strand"], int(r["start"]), int(r["end"])),
        )
    )

    assigned = sum(len(c) for c in clusters)
    if assigned != len(rows):
        sys.exit(f"orfmerge: clustering dropped rows ({assigned} of {len(rows)} assigned)")

    write_catalogue(prefix, clusters, bed_index, min_callers=args.min_callers, min_samples=args.min_samples)
    write_versions()
    return 0


sys.exit(main())
