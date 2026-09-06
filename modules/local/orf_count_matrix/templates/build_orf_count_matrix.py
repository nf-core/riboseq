#!/usr/bin/env python3
"""Pivot per-sample ORF P-site count TSVs into a single ORF x sample matrix.

Per-sample input is `sample_id<TAB>orf_id<TAB>count` (sample_id is identical
on every row, prepended upstream by the per-sample counter); a sample with
zero in-frame P-sites produces an empty file with no sample_id at all, so
the expected sample list is passed in separately rather than inferred from
the concatenated rows. Output rows follow the BED12 catalogue's 4th column
in catalogue order, zero-filled for any ORF or sample with no matching
counts, so the matrix is keyed on the catalogue and the full sample list.
"""

import platform
from pathlib import Path

import pandas as pd
import yaml

counts = pd.concat(
    pd.read_csv(p, sep="\\t", names=["sample", "orf_id", "count"]) for p in sorted(Path("counts").iterdir())
)

catalogue_orfs = (
    pd.read_csv("$orf_catalogue_bed12", sep="\\t", header=None, comment="#", usecols=[3])[3].drop_duplicates().tolist()
)

expected_samples = sorted(s for s in "$sample_ids_csv".split(",") if s)

matrix = (
    counts.pivot_table(index="orf_id", columns="sample", values="count", aggfunc="sum", fill_value=0)
    .reindex(index=catalogue_orfs, columns=expected_samples, fill_value=0)
    .astype(int)
)
matrix.to_csv("${prefix}.tsv", sep="\\t")

with open("versions.yml", "w") as f:
    yaml.safe_dump(
        {"${task.process}": {"python": platform.python_version(), "pandas": pd.__version__}},
        f,
    )
