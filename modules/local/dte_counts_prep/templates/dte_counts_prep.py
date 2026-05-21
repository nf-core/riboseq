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

The two sample sets must be disjoint; the union of column names
becomes the column space of the output, in `primary then secondary`
order.

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
                       columns. First mapping per primary id wins.

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
    df = pd.read_csv(path, sep="\\t")
    for col in (primary_col, secondary_col):
        if col not in df.columns:
            raise SystemExit(f"{path}: expected `{col}` header")
    df = df[[primary_col, secondary_col]].dropna()
    df = df[(df[primary_col] != "") & (df[secondary_col] != "")]
    df = df.drop_duplicates(subset=[primary_col], keep="first")
    return df.set_index(primary_col)[secondary_col]


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
    mapping = load_mapping("${mapping}", opts.primary_id_col, opts.secondary_id_col)

    overlap = set(primary.columns) & set(secondary.columns)
    if overlap:
        raise SystemExit(
            "Primary and secondary matrices share sample columns: "
            + ", ".join(sorted(overlap))
        )

    keep_primary = primary.index[primary.index.isin(mapping.index)]
    mapped_secondary = mapping.loc[keep_primary]
    keep_primary = keep_primary[mapped_secondary.isin(secondary.index)]
    mapped_secondary = mapping.loc[keep_primary]

    primary_kept = primary.loc[keep_primary]
    secondary_expanded = secondary.loc[mapped_secondary.values]
    secondary_expanded.index = keep_primary

    combined = pd.concat([primary_kept, secondary_expanded], axis=1)

    if opts.integer:
        combined = combined.round().astype(int)

    combined.index.name = opts.primary_id_col
    combined.reset_index().to_csv("${prefix}.tsv", sep="\\t", index=False)

    n_dropped_no_map = len(primary.index.difference(mapping.index))
    n_dropped_no_secondary = len(primary.index) - n_dropped_no_map - len(keep_primary)
    sys.stderr.write(
        f"Primary rows: {len(primary.index)}; "
        f"dropped {n_dropped_no_map} with no mapping, "
        f"{n_dropped_no_secondary} mapped but absent from secondary; "
        f"kept {len(keep_primary)}.\\n"
    )

    write_versions()
    return 0


sys.exit(main())
