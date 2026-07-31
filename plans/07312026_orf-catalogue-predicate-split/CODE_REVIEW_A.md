# Code Review A — `feat/orf-catalogue-without-novel-source` (commit `9f4a929`)

Reviewer A. Reviewed against `origin/dev`. All claims below were verified against the code,
the config, the committed snapshots, or a Nextflow probe run — not against the commit message
or the plan.

## Summary

The core redesign is correct. The construction hazard that killed Revision 1 of the plan is
genuinely gone, and I could not construct any configuration in which a `.out` reference or a
`channel.empty()` is reached in a newly-enabled path. The predicate split preserves caller
routing exactly. The dead-code removal is verified safe.

The defect is in the **new test**, which cannot pass as written and has never been run:
one negative assertion is false by construction, and no snapshot file is committed.

| Priority | Count | Area |
|---|---|---|
| Critical | 1 | new test asserts a directory is absent that Nextflow always creates |
| High | 1 | no committed `.snap` for the new test |
| Medium | 3 | denominator namespace, predicate triplication, matrix-shape divergence |
| Low | 6 | redundant cast, stale comments, changelog ordering, overstated claim |

One Critical fix applied. Everything else is a recommendation.

---

## Verified correct (no action)

Recorded so the next reader does not re-derive it.

**No construction hazard in any reachable configuration.** I enumerated the four predicate
combinations and audited every `.out` reference and every channel that can be
`channel.empty()`:

- `workflows/riboseq/main.nf:431` — `if (extended_orf_active || orf_catalogue_active)` is a
  superset of every gate that reads `BUILD_FULL_HYBRID_GTF.out` (`:692`, `:699`, both nested
  inside `if (extended_orf_active)` at `:690`) or `ch_full_hybrid_gtf` (`:472`, `:707`).
- `workflows/riboseq/main.nf:462`, `:531`, `:684` all gate on the same `orf_catalogue_active`,
  so `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.*` (`:536`, `:754`, `:780`) is never read
  un-invoked.
- `ch_orf_count_matrix` (`:506`) is assigned at `:540` under `!skip_plastid && orf_catalogue_active`
  and consumed at `:752` under exactly the same pair.
- `QUANTIFY_STAR_SALMON` is invoked at `:551`, before the new `else` branch reads its output at
  `:748`. Construction order is fine.

**Newly populating `ch_full_hybrid_gtf` does not leak into the callers.**
`subworkflows/local/orf_caller_dispatch/main.nf:63-65` puts `ch_full_hybrid_gtf` behind an
`extended_orf_active` ternary, so in the new configuration (`orf_catalogue_active` true,
`extended_orf_active` false) Rp-Bp (`:166`) and PRICE (`:184`) still receive `ch_gtf`, as
before. This was the main behaviour-change risk of widening the guard at `:431`, and it does
not materialise. The updated contract comment at `:43` is accurate.

**`ch_hybrid_gtf` is always a usable value channel at `:433`.** `novel_source_configured`
(`:370`) is algebraically identical to the novel-discovery gate at `:324`
(`(!skip_stringtie && !novel_gtf) || novel_gtf` ≡ `!skip_stringtie || novel_gtf`), so when
`extended_orf_active` is false the discovery block provably did not run and `ch_hybrid_gtf` is
still `ch_canonical_gtf` from `:311`. No divergence to exploit.

**Superset invariant holds.** `assets/build_full_hybrid_gtf.awk:13-26` prints every row of
file 1 (`ch_gtf`), then only file-2 rows whose `transcript_id` (`:24`) or, failing that,
`gene_id` (`:26`) is absent. With no novel source, file 2 is `ch_canonical_gtf`, so the output
is `ch_gtf` plus canonical-only transcripts. The plan's claim is correct.

**Dead code removal is safe.** `git grep` on `origin/dev` finds
`orf_agreement_min_callers`, `rank_aggregation_callers`, `ch_enabled_orf_callers` and
`ch_rank_aggregation_callers` only at their four definition sites
(`origin/dev:workflows/riboseq/main.nf:460,463,466,467`). No consumers.

**`dotseqPrerequisitesError` cannot emit an empty requirements list.** The gate at
`utils_nfcore_riboseq_pipeline/main.nf:257` fires only if `!orf_catalogue_active ||
skip_plastid || !contrasts`; `!orf_catalogue_active` implies `!extended_orf_analysis ||
!any_caller_enabled`, each of which appends to `missing` (`:259`, `:260`). The empty-list
regression the plan warned about is avoided.

**Plastid warning now fires.** Dropping the `novel_source_configured` conjunct at `:199` makes
`--extended_orf_analysis true --skip_stringtie true --skip_plastid true` warn, which
`docs/usage.md` promises. `.count(true)` at `:207` is correct — every element is a `Boolean`.

**Re-indented block is semantically correct.** `ch_orf_dte_rna_counts` is assigned *without*
`def` in both branches (`:744`, `:748`), so it lands in workflow binding scope and is visible
at `:753`. A `def` there would have been a block-scoped variable and an immediate
`MissingProperty` at `:753`. The code gets this right.

**`nextflow lint` is clean.** The only error in the whole repo is pre-existing nf-core
template code (`subworkflows/nf-core/utils_nextflow_pipeline/main.nf:77`,
`groovy.json.JsonGenerator.Options`). No new warnings attributable to this change.

**Asserted output paths are right.** `orf_catalogue/cohort.catalogue.bed12` matches
`conf/modules.config:1104-1106` and `tests/dotseq.nf.test.snap:819`;
`orf_quantification/orf_psite_counts.tsv` matches `conf/modules.config:1131` and
`tests/dotseq.nf.test.snap:867`. `!hybrid_dir.exists()` is correct — every `hybrid_star`
publisher (`conf/modules.config:1363-1460`) is under `EXTENDED_ORF_SECOND_PASS_ALIGN`, which
is not invoked when `extended_orf_active` is false.

---

## Issues

### Critical

#### C1 — `assert !salmon_hybrid.exists()` is false by construction (FIXED)

`tests/orf_catalogue_no_novel_source.nf.test:33,43` (as committed):

```groovy
def salmon_hybrid = file("${params.outdir}/quantification/salmon_hybrid")
...
{ assert !salmon_hybrid.exists() },
```

**Nextflow creates a `publishDir` target directory whenever the process runs, even if `saveAs`
returns `null` for every output file.** `BUILD_FULL_HYBRID_GTF` now runs in exactly this
configuration — that is the point of the change — and its publish path is
`${params.outdir}/quantification/salmon_hybrid` (`conf/modules.config:716`) with `saveAs`
returning `null` unless `--save_reference` (`:718`). So the directory exists, empty, and the
assertion fails. The test cannot pass.

Three independent confirmations:

1. Probe run on the local Nextflow 25.10.4 — a process whose `saveAs` returns `null` for its
   only output still leaves the publish directory behind:
   `results_probe/never_published` exists and is empty.
2. `tests/default.nf.test.snap:366-368` lists `"genome"`, `"genome/index"` and
   `"genome/sortmerna"` with **no files beneath any of them**. Every publisher to those paths
   (`conf/modules.config:28-146`) is `--save_reference`-gated, and the default test does not
   set it. `genome/index_te` is absent because `KALLISTO_INDEX_TE` did not run — so the
   directory tracks *process execution*, not publication.
3. `tests/dotseq.nf.test.snap:1621` lists `"quantification/salmon_hybrid/star_index"` as an
   empty directory, from the same `--save_reference`-gated pattern at
   `conf/modules.config:734-736`.

This also means the plan's Context bullet "`BUILD_FULL_HYBRID_GTF` publishes nothing unless
`--save_reference` … so running it more often has no snapshot impact" is wrong in a way that
matters: it publishes no *files*, but it does add an empty directory to `stable_name`. Harmless
for the existing snapshots (all extended tests already run it), but it is what breaks C1.

**Fix applied** — assert on the products the hybrid path actually creates rather than on the
shared directory:

```groovy
// BUILD_FULL_HYBRID_GTF creates quantification/salmon_hybrid/ even though it
// publishes nothing without --save_reference, so assert on its contents.
def hybrid_rna_counts = file("${params.outdir}/quantification/salmon_hybrid/salmon.merged.hybrid.gene_counts_length_scaled.tsv")
def hybrid_star_index = file("${params.outdir}/quantification/salmon_hybrid/star_index")
...
{ assert !hybrid_rna_counts.exists() },
{ assert !hybrid_star_index.exists() },
```

Both targets are stronger than the original intent: the first is the denominator matrix itself
(`conf/modules.config:777-783`, present as
`quantification/salmon_hybrid/salmon.merged.hybrid.gene_counts_length_scaled.tsv` at
`tests/dotseq.nf.test.snap:1611`), the second proves no second STAR index was built
(`conf/modules.config:732-737`).

### High

#### H1 — No committed snapshot for the new test

`tests/orf_catalogue_no_novel_source.nf.test:44-51` asserts `snapshot(...).match()`, but there
is no `tests/orf_catalogue_no_novel_source.nf.test.snap`. Every other test in `tests/` has one
(the only exception is the `.disabled` file). CI runs `nf-test test … --ci`
(`.github/actions/nf-test/action.yml:66`, nf-test 0.9.4 per `.github/workflows/nf-test.yml:21`);
`--ci` does not write snapshots, so the assertion either fails outright or is vacuous. Either
way the snapshot half of the test is not doing its job.

This corroborates C1: the commit message says the change was "Verified under `-preview`", which
means the pipeline was never actually executed in this configuration, so neither the snapshot
nor the false assertion could have been caught.

I cannot generate the snapshot here — it needs a full pipeline run with containers and network.
**Run the test once (`nf-test test tests/orf_catalogue_no_novel_source.nf.test`) and commit the
`.snap`.** Expect `quantification/salmon_hybrid` to appear in `stable_name` as an empty
directory; that is correct and should be accepted.

### Medium

#### M1 — The new `else` branch reintroduces a namespace mismatch under `--canonical_gtf`

`workflows/riboseq/main.nf:746-748` feeds `DTE_COUNTS_PREP` with
`QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled`, whose gene ids come from `ch_gtf`
(`:556`). The ORF-to-gene map it is joined against (`:754`) comes from the catalogue, which is
built against `ch_full_hybrid_gtf` (`:472`) — and that is `ch_gtf` **plus any
`--canonical_gtf`-only transcripts** (`assets/build_full_hybrid_gtf.awk:24-26`).

`nextflow.config:27` exposes `--canonical_gtf`, and the schema recommends MANE Select for
human, whose ids are absent from an Ensembl `--gtf`. Those ORFs get a `gene_id` the denominator
matrix has no row for, and
`modules/local/dte_counts_prep/templates/dte_counts_prep.py` drops them silently — the only
signal is the `mapped but no host gene in secondary` count on stderr.

This is the same class of hazard the plan rejected the `ch_gtf` fallback over ("`orfnormalise`
degrades silently"), just moved from the annotation to the denominator. The hybrid branch does
not have it, because `QUANTIFY_HYBRID_RNA` is quantified against `ch_full_hybrid_gtf` (`:733`).

The commit message's justification — "Without one there are no novel genes, so the primary
matrix already covers every host gene" — holds only when `--canonical_gtf` is derived from
`--gtf`, which is the default but not the documented recommendation.

The concurrent working tree already documents this in `docs/usage.md` (see Working-tree note
below), which is the right minimum. **Recommend additionally** a warning when
`params.canonical_gtf` is user-supplied and the `else` branch is taken, or promoting
`dte_counts_prep.py`'s drop count to a `log.warn`-visible failure when it exceeds some fraction
of rows. Docs alone leave a silent-wrong-answer path.

#### M2 — The caller-enabled predicate now exists in three places

- `workflows/riboseq/main.nf:361-366` — `enabled_orf_callers` (a list)
- `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf:204-207` — `enabled_caller_count`
- `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf:254` — `any_caller_enabled`

Three encodings of one fact, in two files, all reading the same five params. The plan
acknowledges the `main.nf`/validator split is hard to unify (validators run before the
workflow), but the two copies *within* the same validator file are trivially unifiable — hoist
one `enabled_caller_count` above both functions, or have `dotseqPrerequisitesError` take it as
an argument. Adding a caller flag currently means editing three sites, and missing one fails
silently rather than loudly.

#### M3 — The two ORF-DTE branches produce different matrix shapes

`dte_counts_prep.py` drops only the secondary columns that overlap the primary (the Ribo-seq
ones). In the hybrid branch the secondary matrix contains RNA-seq columns only, because
`FASTQ_ALIGN_STAR_FULL_HYBRID` is fed `sample_type == 'rnaseq'` reads (`main.nf:702`). In the
new `else` branch the secondary matrix is the whole-samplesheet Salmon matrix, so **TI-seq
columns survive into the combined ORF matrix** alongside the RNA-seq ones.

This is consistent with the gene-level path (`ch_te_counts` at `:573` also carries every sample
type), so it is unlikely to break anything, and the test profile may well have no TI-seq
samples. But the two ORF-DTE branches are no longer shape-equivalent, and the code comment at
`:747` ("`DTE_COUNTS_PREP` drops the Ribo-seq columns it carries") describes only half of what
survives. **Recommend** either filtering the reused matrix to RNA-seq columns for parity, or
extending the comment to say TI-seq columns pass through.

### Low

- **L1** `main.nf:372` — `as boolean` is redundant. Groovy's `&&` already returns a `Boolean`;
  verified by probe (`true && ['ribotish']` → `java.lang.Boolean` / `true`,
  `true && []` → `Boolean` / `false`). The plan's justification ("evaluates to a `List`, which
  is … wrong if it is ever passed as a subworkflow `val`") is factually incorrect. Harmless as
  defensive code — no fix needed, but do not propagate the reasoning.
- **L2** `main.nf:361-366` — after the dead-code deletion, `enabled_orf_callers` is a list of
  strings built solely to be tested for emptiness at `:372`. A count would say the same thing
  and match the validator's shape (M2).
- **L3** `main.nf:859` — emit comment still reads "empty unless extended-ORF + plastid both
  active". The gate is now `orf_catalogue_active && !skip_plastid`. Not fixed here
  deliberately: `main.nf` is currently untouched in the working tree and a formatter hook
  rewrote a whole neighbouring file during this review (see below), so I would not risk a
  900-line reformat for a comment.
- **L4** `tests/dotseq.nf.test:20-24` — "Default ribotish + ribocode is enough" is wrong;
  `conf/test.config:33` sets `skip_ribocode = true`, so the test is single-caller. Plan item 11
  asked for this; not done.
- **L5** `CHANGELOG.md:46` — `#222` is inserted between `#204` and `#210`, while the tail of the
  section runs `#195, #199, #204, #210, #220`. Move it after `#220`.
- **L6** The commit message's "This also makes the ORF-level denominator exactly the gene-level
  one" overstates. The gene-level matrix is canonical-filtered by `FILTER_COUNTS_CANONICAL`
  (`main.nf:649-654`), and under the default `te_quantification_method = 'plastid_psite'` its
  Ribo-seq columns are replaced by plastid P-site counts (`main.nf:625-632`). The RNA columns
  are the same numbers; the matrices are not the same matrix.

---

## Fixes applied

One file, `tests/orf_catalogue_no_novel_source.nf.test` — C1 above. Left uncommitted. Nothing
else was touched by me.

## Working-tree note (not my change)

`git status` shows two files I did not edit:

- `docs/usage.md` — the stale RNA-denominator paragraph at `:439` (plan item 10, dropped in the
  commit) has been rewritten, and the `--canonical_gtf` namespace hazard of M1 documented at
  `:658`. Both are correct and I have not duplicated them. One nit: the new text repeats the
  L6 overstatement ("The ORF-level denominator is then exactly the gene-level one").
- `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf` — **175 lines changed**. A
  formatter has rewritten the entire file: include-block alignment, `take:` comment spacing,
  `UTILS_NEXTFLOW_PIPELINE (` → `(`, trailing commas, `log.warn "…"` → `log.warn("…")`,
  single-line `if` bodies expanded to braces. Buried in it is one intentional addition — a
  warning when `extended_orf_analysis` is set with zero callers enabled (a real gap I had also
  flagged; the catalogue is silently skipped today). **This must be reduced to a minimal diff
  before committing.** A whole-file reformat in a behavioural PR is exactly what the repo's
  contribution guidance says not to ship, and it will bury the real change from reviewers.
  `nextflow lint` passes on the reformatted file, so the reformat is at least not breaking.

## Recommendations, in order

1. **H1** — run the new test and commit its `.snap`. Nothing else validates this change end to
   end, and C1 shows what slips through when the test is never executed.
2. **Working-tree note** — revert the whole-file reformat of the validator, keeping only the
   zero-callers warning.
3. **M1** — add a code-level signal (warning, or a loud drop count) for the
   `--canonical_gtf` / primary-matrix namespace mismatch. Docs are necessary but not
   sufficient for a silent-drop path.
4. **M3** — decide whether TI-seq columns should reach the ORF-level matrix, and either filter
   or document.
5. **M2** — unify the two caller-count encodings inside the validator file.
6. **L3, L4, L5** — one-line cleanups, batch them.
