# Plan Coverage Report

**Mode:** B
**Plan(s):** PLAN.md (rev 2)
**Date:** 2026-07-31
**Verdict:** INCOMPLETE — 14 items unresolved (10 MISSING, 4 PARTIAL), spanning 9 distinct gaps

Audited against branch `feat/orf-catalogue-without-novel-source`, single commit `9f4a929`, plus the
current working tree (clean apart from untracked `plans/`).

## Summary

- Total items: 52
- DONE: 38
- PARTIAL: 4
- MISSING: 10
- DEVIATED: 0

Breakdown by ledger section:

| Section | Items | DONE | PARTIAL | MISSING |
|---|---|---|---|---|
| A. Implementation outline (11 steps, expanded to sub-items) | 28 | 23 | 0 | 5 |
| B. Behavior requirements + edge cases | 14 | 13 | 0 | 1 |
| C. Validation checks | 10 | 2 | 4 | 4 |

Every **code-level** step of the outline (steps 1–9) is fully implemented. All gaps sit in
step 10 (docs/comments) sub-items, step 11 (test) sub-items, and validation execution.

Note: item A18 and item B14 are the *same* underlying gap (the cross-caller annotation
asymmetry note), counted once in each section because the plan lists it as both a step-10
sub-item and a Behavior edge case. Distinct gaps: 9.

## Coverage ledger

### A. Implementation outline

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| A1 | Hoist `enabled_orf_callers` above the predicates; define all three predicates in one block; `orf_catalogue_active` cast `as boolean` | Step 1 | DONE | `workflows/riboseq/main.nf:361-372`. Caller list at `:361-366`, predicates at `:370-372`, cast present at `:372` |
| A2 | Delete dead code `rank_aggregation_callers`, `orf_agreement_min_callers`, `ch_enabled_orf_callers`, `ch_rank_aggregation_callers` | Step 2 | DONE | Whole block incl. its comment header removed. Repo-wide grep over `*.nf`/`*.config`/`*.json`/`*.py`/`*.md` returns hits only inside `plans/` |
| A3 | Widen hybrid-GTF guard to `if (extended_orf_active \|\| orf_catalogue_active)` | Step 3 | DONE | `main.nf:431`; body unchanged; `ch_full_hybrid_gtf` assigned at `:437` |
| A4 | Catalogue gate (old `:477`) → `if (orf_catalogue_active)` | Step 4 | DONE | `main.nf:462` |
| A5 | ORF P-site quantification gate (old `:546`) → `if (orf_catalogue_active)`, still nested in the plastid block | Step 4 | DONE | `main.nf:531`, nested inside `if (!params.skip_plastid)` at `:508` |
| A6 | ORF-DTE gate (old `:699`) → `if (orf_catalogue_active && !params.skip_plastid)` | Step 4 | DONE | `main.nf:684` |
| A7 | No change at old `:487` — catalogue still receives `ch_full_hybrid_gtf` | Step 4 | DONE | `main.nf:472` passes `ch_full_hybrid_gtf` to `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE` |
| A8 | Without a novel source, skip `GFFREAD_FULL_HYBRID`, `STAR_GENOMEGENERATE_FULL_HYBRID`, `FASTQ_ALIGN_STAR_FULL_HYBRID`, `BAM_DEDUP_UMI_HYBRID_RNA`, `QUANTIFY_HYBRID_RNA` | Step 5 | DONE | All five sit inside `if (extended_orf_active)` at `main.nf:690-745`. `-preview` log for `--extended_orf_analysis true --skip_stringtie true` contains zero mentions of any of them |
| A9 | Feed `DTE_COUNTS_PREP` with `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled` in the fallback | Step 5 | DONE | `main.nf:747-748` sets `ch_orf_dte_rna_counts`; consumed at `:751-756`. `QUANTIFY_STAR_SALMON` is invoked unconditionally at `:551`, so `.out` is always defined |
| A10 | Fix the `dotseqPrerequisitesError` gate: drop the two predicate defs, add `orf_catalogue_active`, gate on `!orf_catalogue_active \|\| skip_plastid \|\| !contrasts`, make old `:254` a plain `if`, remove unused `novel_source_configured` | Step 6 | DONE | `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf:254-266`. All five sub-parts present; the `else if` is gone, so `missing` can no longer be empty when the gate fires |
| A11 | Drop the `novel_source_configured` conjunct from the plastid warning guard | Step 7 | DONE | `utils_nfcore_riboseq_pipeline/main.nf:199`; wording at `:200` unchanged as the plan required |
| A12 | Reword the `--extended_orf_analysis`-without-novel-source warning; must not claim the catalogue is skipped | Step 8 | DONE | `utils_nfcore_riboseq_pipeline/main.nf:196-198`. New text scopes the no-op to caller annotation and states "the ORF catalogue is still built" |
| A13 | Warn when `params.orf_min_callers` exceeds the enabled caller count | Step 9 | DONE | `utils_nfcore_riboseq_pipeline/main.nf:203-210`. Guarded additionally by `params.extended_orf_analysis` and `enabled_caller_count > 0`; both narrow the warning to reachable configurations rather than changing its content |
| A14 | `nextflow_schema.json` — fix `extended_orf_analysis` `description` and `help_text` | Step 10 | DONE | Both updated; the "Requires --skip_stringtie false or --novel_gtf" sentence is gone and the catalogue is described as not needing a novel source |
| A15 | `docs/usage.md:625` — the "no-op" claim | Step 10 | DONE | Rewritten at `docs/usage.md:625` |
| A16 | `docs/usage.md:439` — RNA denominator provenance in the fallback | Step 10 | **MISSING** | Paragraph unchanged. See Gap 1 |
| A17 | Add notes on coordinate-only comparability and the single-caller consensus caveat | Step 10 | DONE | `docs/usage.md:627` ("Two caveats when comparing runs…") covers both |
| A18 | Add a note on the cross-caller annotation asymmetry | Step 10 | **MISSING** | Nothing added. See Gap 2 |
| A19 | `docs/output.md:493`, `docs/usage.md:645`, `:654` — no edit needed | Step 10 | DONE | Confirmed untouched by `9f4a929` |
| A20 | Stale comment `subworkflows/local/quantify_orf_psite/main.nf:12` | Step 10 | DONE | Now "`--extended_orf_analysis = true` AND a caller ran AND plastid is enabled" |
| A21 | Stale comment `subworkflows/local/orf_caller_dispatch/main.nf:43-44` | Step 10 | DONE | "empty and unread otherwise" replaced with "populated whenever the hybrid GTF is built" |
| A22 | Stale comment `workflows/riboseq/main.nf:870` (`orf_count_matrix` emit) | Step 10 | **MISSING** | Unchanged, now at `main.nf:859`. See Gap 3 |
| A23 | `CHANGELOG.md` — one `[#NNN]` bullet under `## v1.3.0dev` | Step 10 | DONE | `[#222]` bullet added under `## v1.3.0dev` → `### Fixed` (line 46) |
| A24 | New test: `extended_orf_analysis = true`, `skip_stringtie = true`, no `novel_gtf`, keeping `-profile test`'s `contrasts` | Step 11 | DONE | `tests/orf_catalogue_no_novel_source.nf.test`; `novel_gtf` deliberately unset and commented as such; `contrasts` not overridden |
| A25 | Follow the `stable_name`/`stable_path` + `tests/.nftignore_orf` pattern from `tests/novel_gtf.nf.test:28-46` | Step 11 | DONE | Identical shape; `tests/.nftignore_orf` exists and is referenced |
| A26 | Assert non-zero sizes, not mere existence | Step 11 | DONE | `catalogue_bed12.size() > 0` and `orf_counts.size() > 0`. Asserted paths match the names real runs publish (`orf_catalogue/cohort.catalogue.bed12`, `orf_quantification/orf_psite_counts.tsv`, both present in `tests/stringtie_extended.nf.test.snap`) |
| A27 | Fix the stale comment at `tests/dotseq.nf.test:20-24` ("Default ribotish + ribocode is enough") | Step 11 | **MISSING** | `tests/dotseq.nf.test` is not touched by `9f4a929`; the comment still reads "Default ribotish + ribocode is enough". See Gap 4 |
| A28 | New test's snapshot | Step 11 | **MISSING** | No `tests/orf_catalogue_no_novel_source.nf.test.snap`; every other pipeline test in `tests/` has one. See Gap 5 |

### B. Behavior requirements and edge cases

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| B1 | Prerequisite is `extended_orf_analysis` + ≥1 caller, with no requirement on `skip_stringtie`/`novel_gtf` | Behavior 1 | DONE | `main.nf:372` |
| B2 | Catalogue always built against `ch_full_hybrid_gtf`; no second GTF channel, no fallback ternary | Behavior 2 | DONE | Only `ch_full_hybrid_gtf` exists; no `ch_catalogue_gtf` anywhere |
| B3 | `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE` runs, emitting bed12 / tsv / ORF-to-gene / AA FASTA as today | Behavior 3 | DONE | Process invocation at `main.nf:469-474` unchanged apart from the enclosing gate |
| B4 | `QUANTIFY_ORF_PSITE` runs when the catalogue exists and `--skip_plastid` is false | Behavior 4 | DONE | `main.nf:531` inside `:508` |
| B5 | ORF-level DTE runs when catalogue exists, plastid on, `--contrasts` supplied | Behavior 5 | DONE | `main.nf:684` inside `if (ch_contrasts_file)` at `:640` |
| B6 | Caller routing untouched without a novel source — no second STAR pass, no hybrid transcriptome, no hybrid GTF to Ribo-TISH | Behavior 6 | DONE | `main.nf:376` still `if (extended_orf_active)`; `extended_orf_active` still the dispatch arg at `:449` |
| B7 | Superset invariant holds by construction | Behavior 7 | DONE | Catalogue annotation is `BUILD_FULL_HYBRID_GTF` output; `assets/build_full_hybrid_gtf.awk` unmodified |
| B8 | ORF IDs are not stable across runs; comparability is coordinate-only | Behavior 8 | DONE | Documented at `docs/usage.md:627` |
| B9 | Single caller — catalogue builds; warn rather than document | Edge case | DONE | Warning at `utils_nfcore_riboseq_pipeline/main.nf:208-210`; fired in the `--orf_min_callers 2` preview |
| B10 | No callers enabled — no catalogue, unchanged | Edge case | DONE | `orf_catalogue_active` is `false` for an empty caller list via the `as boolean` cast at `:372` |
| B11 | `--extended_orf_analysis false` — no catalogue; default runs unchanged | Edge case | DONE | Both predicates carry `params.extended_orf_analysis` |
| B12 | `--skip_plastid true` — catalogue builds, no `orf_quantification/`, no ORF-DTE, and the warning fires | Edge case | DONE | Guard fixed at `:199`; warning observed in `.nextflow.log.1` (`--extended_orf_analysis true --skip_stringtie true --skip_plastid true`) |
| B13 | Ribo-seq-only samplesheet + `--contrasts` becomes a clear error via `dte_counts_prep` | Edge case | DONE | Follows from A9: `modules/local/dte_counts_prep/templates/dte_counts_prep.py:130-143` raises `SystemExit` when the secondary matrix has no sample columns left after dropping primary-role overlaps. No code change was required, and no test exercises it |
| B14 | Cross-caller annotation asymmetry — docs only | Edge case | **MISSING** | Same gap as A18. See Gap 2 |

### C. Validation

Evidence comes from the six `.nextflow.log*` files dated 31 Jul 11:35–11:40, all of which are
`nextflow run . -profile test,docker … -preview` invocations. **No full pipeline execution and no
`nf-test` run is recorded**, so every filesystem-level assertion in the plan is unverified.

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| C1 | Changed path works: `--extended_orf_analysis true --skip_stringtie true`, plastid on — hybrid GTF built, catalogue runs, `orf_quantification/` created, no `EXTENDED_ORF_SECOND_PASS_ALIGN`, non-empty bed12/tsv | Validation 1 | PARTIAL | `-preview` (`.nextflow.log.5`) constructs cleanly and mentions `BUILD_FULL_HYBRID_GTF`, `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE`, `QUANTIFY_ORF_PSITE`, `DTE_COUNTS_PREP` with zero `EXTENDED_ORF_SECOND_PASS` mentions. Construction only — nothing executed, so the non-empty-output assertions are unproven |
| C2 | No regression with a novel source: `tests/stringtie_extended.nf.test` passes with an unchanged snapshot | Validation 2 | PARTIAL | Only `-preview` with `--skip_stringtie false --extended_orf_analysis true` (`.nextflow.log.3`), which constructs and shows the hybrid quant processes still present. The nf-test and its snapshot comparison were not run |
| C3 | Fallback annotation equals the reference: `--save_reference` run, diff `BUILD_FULL_HYBRID_GTF` output against `--gtf` | Validation 3 | **MISSING** | No `--save_reference` run in any log; no diff performed |
| C4 | Catalogue comparability on coordinates: two runs differing only by `--novel_gtf`, compare on `chrom/start/end/strand` | Validation 4 | **MISSING** | Requires two full runs; none exist |
| C5 | `--canonical_gtf` superset holds: non-empty `gene_id`, multi-block BED12 for multi-exon ORFs | Validation 5 | **MISSING** | No `--canonical_gtf` run. As the plan's own risk (c) anticipated, the superset invariant currently rests on reading the awk |
| C6 | dotseq no longer over-validates — proceeds, and specifically no empty-requirements error | Validation 6 | DONE | `.nextflow.log.2` (`--extended_orf_analysis true --skip_stringtie true --translational_efficiency_method dotseq`) completes construction with no error. `dotseqPrerequisitesError()` runs during parameter validation, which `-preview` does execute, so this check is genuinely exercised. The remainder of the run was not |
| C7 | Default runs unchanged: `tests/default.nf.test` passes with an unchanged snapshot, no catalogue | Validation 7 | PARTIAL | `-preview` of plain `-profile test,docker` (`.nextflow.log.4`) constructs with no warnings. The nf-test and snapshot comparison were not run |
| C8 | Warnings fire: reworded `:194` warning, plastid warning, new `orf_min_callers` warning | Validation 8 | DONE | All three observed. Reworded warning in `.nextflow.log.5`/`.1`/`.2`; plastid warning in `.nextflow.log.1`; `--orf_min_callers 2` warning in `.nextflow.log` |
| C9 | No redundant hybrid quant: `quantification/salmon_hybrid/` absent in the new test | Validation 9 | PARTIAL | Construction-level only: zero `QUANTIFY_HYBRID_RNA`/`STAR_GENOMEGENERATE_FULL_HYBRID`/`GFFREAD_FULL_HYBRID` in `.nextflow.log.5` versus 20/5/5 mentions in the novel-source preview `.nextflow.log.3`. The test's `!salmon_hybrid.exists()` assertion has never been run |
| C10 | Linting: `nf-core lint` (schema text changed) and `nextflow lint` (removed variables, new predicate) | Validation 10 | **MISSING** | No evidence of either being run. Related observation: the `-preview` runs that build the hybrid GTF emit `WARN … The operator 'first' is useless when applied to a value channel … check channel 'ch_gtf'`. It also appears in the novel-source preview, so it is pre-existing on this path, but it now surfaces in catalogue-only runs too |

## Gaps (detail)

### Gap 1 — `docs/usage.md:439` still describes the pre-change RNA denominator (item A16)

**Expected:** Step 10 required `docs/usage.md:439` to be updated because "RNA denominator provenance
changes in the fallback (relevant if step 5 lands)". Step 5 landed.

**Found:** The paragraph is unchanged. It states that the ORF-level RNA-seq matrix "is quantified
against the full reference transcriptome augmented with the novel intergenic transcripts, not the
canonical reference alone, using the same STAR-alignment then Salmon path as the primary RNA-seq
quantification", and that "The augmented RNA-seq matrices are published under
`<outdir>/quantification/salmon_hybrid/`".

**Gap:** Both sentences are now false whenever there is no novel source — the fallback reuses the
primary Salmon matrix (`main.nf:748`) and publishes no `salmon_hybrid/` directory, which the new
test asserts. The paragraph needs a fallback clause.

### Gap 2 — cross-caller annotation asymmetry is not documented (items A18, B14)

**Expected:** Step 10 required a note on the cross-caller annotation asymmetry, and the Behavior
section marks the same case "Docs only": in fallback mode Ribo-TISH and Ribotricer read the
one-transcript-per-gene canonical GTF while Rp-Bp, PRICE and RiboCode read the multi-isoform
`ch_gtf`, so an isoform-specific ORF is invisible to the former and can never reach 2-of-2
agreement.

**Found:** Nothing. The added paragraph at `docs/usage.md:627` covers the other two notes
(comparability, single-caller consensus) but not this one. The nearest existing text,
`docs/usage.md:387` and `:399`, explains why Rp-Bp and PRICE get the full multi-isoform GTF; it
does not state the consequence for cross-caller agreement.

**Gap:** Add the asymmetry note to the extended-ORF section of `docs/usage.md`.

### Gap 3 — stale `orf_count_matrix` emit comment (item A22)

**Expected:** Step 10 listed `workflows/riboseq/main.nf:870` among the stale gating comments to
update.

**Found:** Unchanged, now at `main.nf:859`: `// channel: [ meta, orf_psite_counts.tsv ] - per-ORF
P-site count matrix; empty unless extended-ORF + plastid both active`. Verified that pre-change
line 870 is exactly this line (`git show 9f4a929^:workflows/riboseq/main.nf`).

**Gap:** The comment's gating description no longer matches `main.nf:531` (`orf_catalogue_active`
= flag + ≥1 caller). Update it.

### Gap 4 — stale comment in `tests/dotseq.nf.test` (item A27)

**Expected:** Step 11 explicitly required fixing the stale comment at `tests/dotseq.nf.test:20-24`
("Default ribotish + ribocode is enough" — ribocode is off in the test profile).

**Found:** `tests/dotseq.nf.test` is not in the commit's file list. Lines 18-23 still read
"… Default ribotish + ribocode is enough to supply the ORF catalogue that DOTSeq consumes."

**Gap:** Correct the comment to reflect that `-profile test` sets `skip_ribocode = true`
(`conf/test.config:33`), so the test is single-caller.

### Gap 5 — new test has no recorded snapshot (item A28)

**Expected:** The new test calls `snapshot(removeNextflowVersion(...), stable_name,
stable_path).match()`, following the `novel_gtf` pattern, whose snapshot is committed as
`tests/novel_gtf.nf.test.snap`.

**Found:** No `tests/orf_catalogue_no_novel_source.nf.test.snap`. Every other pipeline test in
`tests/` ships a `.snap`.

**Gap:** Run the test to record the snapshot and commit it. Until then the test has never executed,
so all five of its assertions — including the two non-empty-size checks and the two
absence checks that stand in for Validations 1 and 9 — are unverified, and a CI run in `--ci` mode
would fail on the missing snapshot.

### Gap 6 — Validations 3, 4, 5 and 10 not performed (items C3, C4, C5, C10)

**Expected:** A `--save_reference` diff proving the fallback annotation equals `--gtf`; a two-run
coordinate comparison; a `--canonical_gtf` superset run asserting non-empty `gene_id` and
multi-block BED12; `nf-core lint` and `nextflow lint`.

**Found:** No evidence for any of the four. All six recorded Nextflow invocations are `-preview`.

**Gap:** Each requires an actual run (or a lint invocation). Validation 5 is the one the plan itself
flagged as load-bearing: without it, the superset invariant that protects `orfnormalise.py` from
silently emitting single-block BED12s rests on code reading alone.

## Test verification

| Test | File | Status |
|---|---|---|
| `-profile test --extended_orf_analysis true` with no novel-transcript source | `tests/orf_catalogue_no_novel_source.nf.test` | PRESENT, NEVER RUN — no `.snap`, no `nf-test` invocation recorded |
| `workflow.success` | same | not exercised |
| `orf_catalogue/cohort.catalogue.bed12` exists and non-empty | same | not exercised (path name matches `stringtie_extended` snapshot, so it is well-formed) |
| `orf_quantification/orf_psite_counts.tsv` exists and non-empty | same | not exercised (path name matches) |
| `hybrid_star/` absent | same | not exercised. Path is correct: `conf/modules.config:1363-1460` publishes every `EXTENDED_ORF_SECOND_PASS_ALIGN` output under `${params.outdir}/hybrid_star` |
| `quantification/salmon_hybrid/` absent | same | not exercised. Path is correct: `conf/modules.config:716-788` |
| Snapshot match | same | MISSING — no snapshot file recorded |
| `tests/stringtie_extended.nf.test` (Validation 2) | existing | not run on this branch |
| `tests/default.nf.test` (Validation 7) | existing | not run on this branch |
| `tests/dotseq.nf.test` | existing | not run on this branch; its stale comment is also unfixed (Gap 4) |

## Verdict

**INCOMPLETE — 14 items unresolved (10 MISSING, 4 PARTIAL), 9 distinct gaps.**

All of the plan's code changes are in place: steps 1–9 are DONE without deviation, all eight
Behavior requirements hold, and five of the six edge cases are covered in code. What remains:

1. **`docs/usage.md:439`** (step 10) — the ORF-level RNA denominator paragraph still claims the
   matrix is quantified against the novel-augmented transcriptome and published under
   `quantification/salmon_hybrid/`. Both are false in the no-novel-source path added by step 5.
2. **Cross-caller annotation asymmetry note** (step 10, Behavior edge case) — not written anywhere
   in `docs/`.
3. **`workflows/riboseq/main.nf:859`** (step 10) — the `orf_count_matrix` emit comment still says
   "empty unless extended-ORF + plastid both active"; the plan listed this line explicitly.
4. **`tests/dotseq.nf.test:18-23`** (step 11) — the "Default ribotish + ribocode is enough" comment
   was to be corrected; the file is untouched.
5. **Record the new test's snapshot** and run it (step 11). Until it runs, Validations 1 and 9 rest
   on `-preview` construction evidence only.
6. **Validation 2 and 7** — run `tests/stringtie_extended.nf.test` and `tests/default.nf.test` to
   confirm the snapshots really are unchanged. The plan's Integration section says "Verify rather
   than assume".
7. **Validation 3** — `--save_reference` run, diff `BUILD_FULL_HYBRID_GTF` output against `--gtf`.
8. **Validation 4** — two runs differing only by `--novel_gtf`, compared on coordinates.
9. **Validation 5** — the `--canonical_gtf` superset run, which the plan identified as the check
   that makes the chosen design provably safe. **Validation 10** — `nf-core lint` and
   `nextflow lint`.
