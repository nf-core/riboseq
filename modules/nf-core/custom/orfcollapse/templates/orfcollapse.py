#!/usr/bin/env python3
"""Collapse small ORFs sharing an amino-acid cluster into a single catalogue entry.

The coordinate-based merge in custom/orfmerge groups ORFs that overlap on the
genome, but the same micropeptide is frequently encoded at several distinct,
non-overlapping genomic loci (typically repetitive regions), and those copies
survive as separate catalogue rows. Following the GENCODE Ribo-seq ORF catalogue
convention (Mudge et al. 2022, Nat Biotechnol, doi:10.1038/s41587-022-01369-0;
gencode-riboseqORFs collapse_cutoff 0.9), small ORFs (orf_class == "smORF", i.e.
aa_length <= 100) are clustered by amino-acid sequence identity upstream
(mmseqs/easycluster) and this module folds each multi-member cluster down to one
representative.

Only smORF rows are collapsed; larger ORFs and transcript-anchored classes pass
through untouched, preserving the deterministic coordinate/transcript merge from
upstream. Among the smORF members of a cluster the representative is chosen here
(longest aa_length, ties broken by orf_id) so the result is independent of which
sequence MMseqs2 labelled the cluster representative. Catalogue row order is
preserved; dropped members fold their cross-caller / cross-sample evidence and
gene mappings into the survivor.
"""

import argparse
import csv
import platform
import re
import shlex
import sys
from collections import OrderedDict, defaultdict

import yaml

CALLERS = ("ribotish", "ribocode", "ribotricer", "rpbp", "price")
# `bedtools getfasta -nameOnly -s` appends the strand as "(+)"/"(-)" to each
# sequence name, so FASTA headers and MMseqs2 cluster ids arrive as
# "<orf_id>(+)". Strip it to recover the bare orf_id used in the catalogue.
STRAND_SUFFIX = re.compile(r"\\([+-]\\)\$")
SCORE_DIRECTIONS = {
    "ribotish": "min",
    "ribocode": "min",
    "ribotricer": "max",
    "rpbp": "max",
    "price": "min",
}
CLASS_ORDER = ("canonical_cds", "uORF", "dORF", "novel_u", "smORF", "other")
SMORF_CLASS = "smORF"


def read_fasta(path):
    seqs = OrderedDict()
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                header = line[1:].split()[0] if len(line) > 1 else ""
                name = STRAND_SUFFIX.sub("", header)
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_clusters(path):
    """Map member orf_id -> cluster key (representative id) from an mmseqs cluster TSV."""
    cluster_of = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) >= 2:
                member = STRAND_SUFFIX.sub("", fields[1])
                rep = STRAND_SUFFIX.sub("", fields[0])
                cluster_of[member] = rep
    return cluster_of


def best_score(values, direction):
    nums = []
    for v in values:
        try:
            if v not in ("", None):
                nums.append(float(v))
        except ValueError:
            continue
    if not nums:
        return ""
    return f"{(max(nums) if direction == 'max' else min(nums)):.6g}"


def merge_members(members):
    """Fold smORF rows sharing an AA cluster into one representative row dict."""
    rep = sorted(members, key=lambda r: (-int(r.get("aa_length") or 0), r["orf_id"]))[0]
    out = dict(rep)
    for c in CALLERS:
        out[f"called_by_{c}"] = (
            "1" if any(r.get(f"called_by_{c}") == "1" for r in members) else "0"
        )
        out[f"score_{c}"] = best_score(
            [r.get(f"score_{c}", "") for r in members], SCORE_DIRECTIONS[c]
        )
    samples = sorted(
        {s for r in members for s in (r.get("samples") or "").split(",") if s}
    )
    out["n_samples"] = str(len(samples))
    out["samples"] = ",".join(samples)
    return rep["orf_id"], out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(shlex.split("${args}"))

    prefix = "${prefix}"

    with open("${catalogue_tsv}", newline="") as fh:
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\\t"
        )
        header = reader.fieldnames
        rows = list(reader)

    bed_index = {}
    with open("${bed12}") as fh:
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) >= 12:
                bed_index[parts[3]] = line.rstrip("\\n")

    with open("${orf_to_gene_tsv}") as fh:
        o2g_header = fh.readline().rstrip("\\n")
        o2g = [line.rstrip("\\n").split("\\t")[:3] for line in fh if line.strip()]

    aa = read_fasta("${aa_fasta}")
    cluster_of = read_clusters("${cluster_tsv}")

    smorf_rows = [r for r in rows if r.get("orf_class") == SMORF_CLASS]
    clusters = defaultdict(list)
    for r in smorf_rows:
        clusters[cluster_of.get(r["orf_id"], r["orf_id"])].append(r)

    remap, merged_rows, dropped = {}, {}, set()
    for members in clusters.values():
        if len(members) < 2:
            continue
        rep_id, merged = merge_members(members)
        merged_rows[rep_id] = merged
        for m in members:
            remap[m["orf_id"]] = rep_id
            if m["orf_id"] != rep_id:
                dropped.add(m["orf_id"])

    per_class = defaultdict(int)
    with (
        open(f"{prefix}.catalogue.tsv", "w") as th,
        open(f"{prefix}.catalogue.bed12", "w") as bh,
        open(f"{prefix}.catalogue.aa.fasta", "w") as ah,
    ):
        th.write("\\t".join(header) + "\\n")
        for r in rows:
            oid = r["orf_id"]
            if oid in dropped:
                continue
            row = merged_rows.get(oid, r)
            th.write("\\t".join(row.get(col, "") for col in header) + "\\n")
            per_class[row.get("orf_class", "other")] += 1
            if oid in bed_index:
                bh.write(bed_index[oid] + "\\n")
            if oid in aa:
                ah.write(f">{oid}\\n{aa[oid]}\\n")

    seen = set()
    with open(f"{prefix}.orf_to_gene.tsv", "w") as oh:
        oh.write(o2g_header + "\\n")
        for orf_id, gene_id, tx_id in o2g:
            orf_id = remap.get(orf_id, orf_id)
            key = (orf_id, gene_id, tx_id)
            if key in seen:
                continue
            seen.add(key)
            oh.write("\\t".join(key) + "\\n")

    with open(f"{prefix}.catalogue.mqc.tsv", "w") as mh:
        mh.write("# id: orf_catalogue\\n")
        mh.write("# section_name: 'ORF catalogue'\\n")
        mh.write("# description: 'Per-class ORF counts in the merged catalogue.'\\n")
        mh.write("# plot_type: 'table'\\n")
        mh.write("# pconfig:\\n")
        mh.write("#   id: 'orf_catalogue_table'\\n")
        mh.write("#   title: 'ORF catalogue'\\n")
        mh.write("Class\\tCount\\n")
        for cls in CLASS_ORDER:
            mh.write(f"{cls}\\t{per_class.get(cls, 0)}\\n")

    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {"${task.process}": {"python": platform.python_version()}},
            fh,
            default_flow_style=False,
            sort_keys=False,
        )
    return 0


sys.exit(main())
