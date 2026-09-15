#!/usr/bin/env python3
"""Combine two count matrices into one feature-level matrix for joint analysis.

The primary matrix is already at the desired feature resolution (rows
indexed by `--primary-id-col`, columns are samples). The secondary
matrix is at a coarser resolution (rows indexed by `--secondary-id-col`,
columns are a different set of samples). A mapping table relates each
primary id to one secondary id.

For every primary row, the secondary row of its mapped id is
column-bound onto the primary row. The result is a single matrix with
one row per primary id and columns = primary samples + secondary
samples. Primary rows whose mapped secondary id is absent from the
secondary matrix are dropped. Primary rows with no mapping in the
mapping table are also dropped.

Primary's column set is authoritative for the primary role. If any
column in secondary overlaps with primary, the secondary copy is
dropped (it almost always represents the same sample quantified in a
different role, e.g. a Salmon all-sample matrix carrying its Ribo-seq
columns alongside the wanted RNA-seq ones). The output column order
is `primary then surviving-secondary`.

When the secondary matrix is expanded by replication (multiple primary
ids share one secondary id), the secondary columns of those rows are
perfectly correlated. Downstream models that assume row independence
across the secondary block will be over-confident; that caveat lives
with the caller of this module.

Inputs (positional, via the Nextflow template):
  ${primary_counts}    TSV: primary id col + per-sample integer columns.
  ${secondary_counts}  TSV: secondary id col + per-sample columns. May
                       carry an extra label column (e.g. `gene_name`)
                       which is auto-dropped if listed in
                       `--secondary-drop-cols`.
  ${mapping}           TSV with at least the primary and secondary id
                       columns. When a primary id maps to several secondary
                       ids, the one whose secondary-matrix row has the highest
                       total count is used.

`ext.args` (parsed below via argparse):
  --primary-id-col STR        Default: `orf_id`.
  --secondary-id-col STR      Default: `gene_id`.
  --secondary-drop-cols STR   Comma-separated columns to drop from the
                              secondary matrix before joining
                              (e.g. `gene_name`). Default: `gene_name`.
  --integer / --no-integer    Cast all count columns to integer.
                              Default: enabled (DESeq2 requires this).

Output:
  `${prefix}.tsv`            TSV with `<primary-id-col>` then
                              primary-sample columns then
                              secondary-sample columns.
"""
import argparse
import platform
import sys

import pandas as pd
import yaml


def load_matrix(path, id_col, drop_cols=None):
    df = pd.read_csv(path, sep="\\t")
    if id_col not in df.columns:
        raise SystemExit(f"{path}: expected `{id_col}` header, got {list(df.columns)}")
    if drop_cols:
        df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    df = df.set_index(id_col)
    return df


def load_mapping(path, primary_col, secondary_col):
    """Return {primary_id: [candidate secondary_ids]}, first-seen order preserved."""
    df = pd.read_csv(path, sep="\\t")
    for col in (primary_col, secondary_col):
        if col not in df.columns:
            raise SystemExit(f"{path}: expected `{col}` header")
    df = df[[primary_col, secondary_col]].dropna()
    df = df[(df[primary_col] != "") & (df[secondary_col] != "")]
    candidates = {}
    for primary_id, secondary_id in zip(df[primary_col], df[secondary_col]):
        genes = candidates.setdefault(primary_id, [])
        if secondary_id not in genes:
            genes.append(secondary_id)
    return candidates


def write_versions():
    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {
                "${task.process}": {
                    "python": platform.python_version(),
                    "pandas": pd.__version__,
                }
            },
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-id-col", default="orf_id")
    parser.add_argument("--secondary-id-col", default="gene_id")
    parser.add_argument(
        "--secondary-drop-cols",
        default="gene_name",
        help="Comma-separated columns to drop from secondary matrix.",
    )
    parser.add_argument("--integer", dest="integer", action="store_true", default=True)
    parser.add_argument("--no-integer", dest="integer", action="store_false")
    raw_args = "${args}".split() if "${args}".strip() else []
    opts = parser.parse_args(raw_args)

    drop_cols = [c.strip() for c in opts.secondary_drop_cols.split(",") if c.strip()]

    primary = load_matrix("${primary_counts}", opts.primary_id_col)
    secondary = load_matrix("${secondary_counts}", opts.secondary_id_col, drop_cols=drop_cols)
    candidates = load_mapping("${mapping}", opts.primary_id_col, opts.secondary_id_col)

    # Sample IDs may overlap when the secondary matrix is wider than the
    # secondary role demands (e.g. a Salmon gene matrix quantified across the
    # whole sample sheet feeds the RNA-seq role, but its columns still include
    # the Ribo-seq samples). Primary's columns are authoritative for the
    # primary role, so drop those columns from secondary rather than failing.
    overlap = sorted(set(primary.columns) & set(secondary.columns))
    if overlap:
        print(
            "Dropping {0} primary-role sample column(s) from secondary: {1}".format(
                len(overlap), ", ".join(overlap)
            ),
            file=sys.stderr,
        )
        secondary = secondary.drop(columns=overlap)
        if secondary.shape[1] == 0:
            raise SystemExit(
                "Secondary matrix has no sample columns left after dropping primary-role "
                "overlaps. Check that the secondary matrix includes the role-specific samples "
                "(e.g. RNA-seq samples for the RNA-seq denominator)."
            )

    # When a (collapsed) ORF maps to several genes, use the host gene present in
    # the secondary matrix with the highest total count as the denominator: it
    # is the gene most likely to drive the ORF's expression, so the better
    # baseline. This uses the raw total (gene length is unavailable here) and is
    # therefore mildly biased towards longer genes, but beats an arbitrary pick.
    gene_totals = secondary.sum(axis=1)
    resolved = {}
    for orf, genes in candidates.items():
        present = [g for g in genes if g in secondary.index]
        if present:
            resolved[orf] = max(present, key=lambda g: gene_totals[g])

    mapping = pd.Series(resolved, name=opts.secondary_id_col)
    keep_primary = primary.index[primary.index.isin(mapping.index)]
    mapped_secondary = mapping.loc[keep_primary]

    primary_kept = primary.loc[keep_primary]
    secondary_expanded = secondary.loc[mapped_secondary.values]
    secondary_expanded.index = keep_primary

    combined = pd.concat([primary_kept, secondary_expanded], axis=1)

    if opts.integer:
        combined = combined.round().astype(int)

    combined.index.name = opts.primary_id_col
    combined.reset_index().to_csv("${prefix}.tsv", sep="\\t", index=False)

    n_no_map = sum(orf not in candidates for orf in primary.index)
    n_no_secondary = sum(
        orf in candidates and orf not in resolved for orf in primary.index
    )
    sys.stderr.write(
        f"Primary rows: {len(primary.index)}; "
        f"dropped {n_no_map} with no mapping, "
        f"{n_no_secondary} mapped but no host gene in secondary; "
        f"kept {len(keep_primary)}.\\n"
    )

    write_versions()
    return 0


sys.exit(main())
