# Decouple the ORF catalogue from the novel-transcript source

> **Revision 2** (2026-07-31). Rewritten after dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). Revision 1 introduced a `ch_catalogue_gtf` ternary while leaving `BUILD_FULL_HYBRID_GTF` gated on `extended_orf_active`; both reviewers independently found that this **crashes at workflow construction** in the configuration the plan exists to enable. See Revision History at the end.

## Goal

Make the cross-caller ORF catalogue — and the ORF quantification and ORF-level DTE that hang off it — depend on `--extended_orf_analysis` alone, rather than on `--extended_orf_analysis` **and** a configured novel-transcript source.

Today the catalogue does not exist unless StringTie ran or `--novel_gtf` was supplied. Since `--skip_stringtie` defaults to `true`, `--extended_orf_analysis true` on its own produces no catalogue, no `orf_quantification/`, and no ORF-level DTE. After this change, extended-ORF mode always produces a catalogue; novel transcripts merely enlarge the annotation it is built against.

Outcome: a run with `--extended_orf_analysis true` and no novel source produces a catalogue over the full reference annotation, **comparable on genomic coordinates** to a novel-source run. Not comparable by `orf_id` — see Behavior item 8.

This is a **code/docs reconciliation, not new behaviour**: `docs/output.md:493`, `docs/usage.md:645` and `docs/usage.md:654` already describe the post-change behaviour ("Produced only when `--extended_orf_analysis true` is set and at least one ORF caller is enabled"). The code is the outlier.

## Context

All affected logic is in `workflows/riboseq/main.nf`, plus parameter validation in `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf`.

One predicate currently serves two unrelated purposes:

```groovy
def novel_source_configured = !params.skip_stringtie || params.novel_gtf   // :361
def extended_orf_active     = params.extended_orf_analysis && novel_source_configured   // :362
```

| Purpose | Needs a novel source? | Sites |
|---|---|---|
| Route hybrid annotation into the callers (second STAR pass, RiboCode hybrid BAM, Ribo-TISH `-g`) | **Yes** — no hybrid transcriptome exists without one | `:366`, `:421`, `:439` |
| Build the catalogue, quantify ORFs, run ORF-level DTE | **No** — needs only caller output plus a GTF | `:477`, `:546`, `:699` |

The approach is to add a second predicate for the lower group and **build the full hybrid GTF whenever either group is active**, rather than giving the catalogue a separate fallback channel. `ch_hybrid_gtf` is always defined — it is `ch_canonical_gtf` at `:311` and is only reassigned under novel discovery at `:341` — so `BUILD_FULL_HYBRID_GTF` runs correctly with no novel source, emitting `ch_gtf` plus any canonical-only transcripts.

Key facts established during review:

- **`:710` and `:717` read `BUILD_FULL_HYBRID_GTF.out.output` directly**, not via `ch_full_hybrid_gtf`. Leaving that process gated while opening the `:699` block raises `Access to 'BUILD_FULL_HYBRID_GTF.out' is undefined…` at construction time. A grep for `ch_full_hybrid_gtf` does **not** surface these lines. Reproduced on Nextflow 25.10.4 by Reviewer B.
- **`enabled_orf_callers` is defined at `:450`, after the `:421` guard** it must now feed. Hoisting it is a prerequisite of this design, not tidy-up.
- **`assets/build_full_hybrid_gtf.awk:13-26`** prints every row of file 1 (`ch_gtf`) verbatim, then only file-2 rows whose `transcript_id`/`gene_id` is absent. With no novel entries the output equals `ch_gtf`; with a user-supplied `--canonical_gtf` it additionally carries the canonical-only transcripts. This is what makes the superset invariant hold by construction.
- **`main.nf:463-467`** (`orf_agreement_min_callers`, `rank_aggregation_callers`, `ch_enabled_orf_callers`, `ch_rank_aggregation_callers`) is **dead code with no consumers** anywhere in the repo. The live consensus threshold is `params.orf_min_callers`, consumed via `conf/modules.config:1081` and `:1104`, default `1` (`nextflow.config:164`).
- **`BUILD_FULL_HYBRID_GTF` publishes nothing** unless `--save_reference` (`conf/modules.config:718`), so running it more often has no snapshot impact. Its publish path is `quantification/salmon_hybrid` (`:716`), which is a misleading location for a fallback run but only visible under `--save_reference`.

## Behavior

1. **Prerequisite** — `params.extended_orf_analysis` is true and at least one ORF caller is enabled. No requirement on `--skip_stringtie` or `--novel_gtf`.
2. **Annotation** — the catalogue is always built against `ch_full_hybrid_gtf`. With a novel source that is reference + novel genes; without one it is the reference plus any canonical-only transcripts. There is no second GTF channel and no fallback ternary.
3. **Catalogue** — `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE` runs, emitting `catalogue_bed12`, `catalogue_tsv`, the ORF-to-gene map and the AA FASTA, as today.
4. **ORF quantification** — `QUANTIFY_ORF_PSITE` runs when the catalogue exists and `--skip_plastid` is false, publishing `orf_quantification/`. Plastid dependency unchanged.
5. **ORF-level DTE** — runs when the catalogue exists, `--skip_plastid` is false and `--contrasts` is supplied.
6. **Caller routing untouched** — without a novel source, callers keep their existing annotations. No second STAR pass, no hybrid transcriptome, no hybrid GTF into Ribo-TISH.
7. **Superset invariant** — the catalogue annotation always contains every transcript ID any enabled caller was given. Guaranteed by the awk, not by assumption. This is load-bearing: violating it makes `orfnormalise.py` silently emit single-block BED12s (`:557-568`, `:599-616`), drop rows (`:473-491`) or produce empty `gene_id`s (`:715`, `:806-819`), which flow into `BEDTOOLS_GETFASTA -split` as wrong sequence and wrong peptides.
8. **ORF IDs are not stable across runs** — `orfmerge.py:206-234` sorts clusters by coordinate and assigns `orf_{idx+1:08d}`, so inserting novel-gene ORFs renumbers everything downstream. Cohort-global `MMSEQS_EASYCLUSTER` + `CUSTOM_ORFCOLLAPSE` can also reshuffle representatives for shared ORFs. Runs are comparable on coordinates only. Users joining two runs' `dte/orf_level/` tables on `orf_id` get a silent systematic mis-join.

Edge cases:

- **Single caller** — catalogue builds. With the default `--orf_min_callers 1` the consensus view equals the full catalogue. With `--orf_min_callers 2` (which `nextflow_schema.json` actively recommends) the consensus view is published **empty apart from headers, silently**. `-profile test` is single-caller (`conf/test.config:33` sets `skip_ribocode = true`), so this is the shape every test exercises. Warn rather than document.
- **No callers enabled** — no catalogue. Unchanged.
- **`--extended_orf_analysis false`** — no catalogue, as now. This change does **not** make the catalogue available in default runs.
- **`--skip_plastid true`** — catalogue builds, no `orf_quantification/`, no ORF-DTE. The warning for this must fire (see step 7); today its guard suppresses it.
- **Ribo-seq-only samplesheet + `--contrasts`** — the RNA-seq presence check at `:316-322` only fires for `run_stringtie`, so the ORF-DTE block is constructed with no RNA-seq input and produces nothing. Step 5 turns this into a clear error via `dte_counts_prep.py:140-142`.
- **Cross-caller annotation asymmetry** — in fallback mode Ribo-TISH and Ribotricer use the one-transcript-per-gene canonical GTF while Rp-Bp, PRICE and RiboCode use multi-isoform `ch_gtf`. An isoform-specific ORF is therefore invisible to the former and can never reach 2-of-2 agreement. Pre-existing in extended mode too, and deliberately not fixed here, but it becomes the first thing a user hits once the catalogue is generally available. Docs only.

## Implementation outline

1. **Hoist the caller set and define all predicates together.** Move `enabled_orf_callers` (`:450-456`) to just above `:361`. It is pure param logic with no channel dependency. Then define, in one block:

   ```groovy
   def novel_source_configured = !params.skip_stringtie || params.novel_gtf
   def extended_orf_active     = params.extended_orf_analysis && novel_source_configured
   def orf_catalogue_active    = (params.extended_orf_analysis && enabled_orf_callers) as boolean
   ```

   `as boolean` matters: `params.extended_orf_analysis && enabled_orf_callers` evaluates to a `List`, which is fine in an `if` but wrong if it is ever passed as a subworkflow `val` (as `extended_orf_active` is at `:439`).

   **This hoist is a prerequisite**, not cleanup: step 3 needs `orf_catalogue_active` at `:421`, which currently precedes the definition at `:450`.

2. **Delete the dead code** at `:460-467` (`rank_aggregation_callers`, `orf_agreement_min_callers`, `ch_enabled_orf_callers`, `ch_rank_aggregation_callers`) — verified to have no consumers. Optional but it removes a trap: it looks like the live consensus threshold and is not.

3. **Widen the hybrid-GTF guard** at `:421`:

   ```groovy
   if (extended_orf_active || orf_catalogue_active) {
   ```

   Nothing else in that block changes. `ch_full_hybrid_gtf` is then always populated whenever any consumer needs it, which is what removes the crash at `:710`/`:717` and the empty-channel hazard at `:725`/`:751`.

4. **Repoint the three catalogue gates**:
   - `:477` — `if (extended_orf_active && enabled_orf_callers)` → `if (orf_catalogue_active)`
   - `:546` — same substitution (stays nested inside `if (!params.skip_plastid)` at `:523`)
   - `:699` — `if (extended_orf_active && enabled_orf_callers && !params.skip_plastid)` → `if (orf_catalogue_active && !params.skip_plastid)`

   No change at `:487` — the catalogue keeps receiving `ch_full_hybrid_gtf`.

5. **Reuse the primary Salmon matrix as the ORF-DTE denominator when there is no novel source.** *In scope, confirmed 2026-07-31.* Inside `:699`, branch on `extended_orf_active`: when false, skip `GFFREAD_FULL_HYBRID` (`:709`), `STAR_GENOMEGENERATE_FULL_HYBRID` (`:715`), `FASTQ_ALIGN_STAR_FULL_HYBRID` (`:722`) and `BAM_DEDUP_UMI_HYBRID_RNA` (`:734`)/`QUANTIFY_HYBRID_RNA` (`:746`), and feed `DTE_COUNTS_PREP` with `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled` instead.

   Justified by `dte_counts_prep.py:17-21` and `:126-142`, which already handle "a Salmon all-sample matrix carrying its Ribo-seq columns alongside the wanted RNA-seq ones" by dropping overlapping primary-role columns. Both matrices are `counts_gene_length_scaled` from the same subworkflow with gene IDs from `ch_gtf`, so shapes match. The rationale at `:700-708` for the separate hybrid quant is explicitly "novel genes have no row in the primary matrix" — void when there are no novel genes.

   Caveat: the primary matrix derives from `ch_transcript_fasta` rather than a gffread extraction of the same GTF, so the two paths are equivalent at the annotation level but not necessarily bit-identical. Acceptable, since the baseline is "no output at all".

   This step is what keeps the change from being a performance regression: without it, a no-novel-source run with contrasts pays for a full STAR genome index build plus a second alignment and Salmon pass to recompute a matrix it already has, and publishes a duplicate `quantification/salmon_hybrid/`. It is also scientifically preferable — the ORF-level RNA denominator becomes literally the same numbers as the gene-level one, which is what the comment at `:704-708` wants.

6. **Fix the `dotseqPrerequisitesError` *gate*, not its message** (`utils_nfcore_riboseq_pipeline/main.nf:244-260`). Deleting only the `missing <<` clause at `:254` leaves the gate at `:251` firing on `!extended_orf_active`, producing `error("… requires: .")` with an empty list — worse than today. Instead: drop `:247-248`, define `def orf_catalogue_active = params.extended_orf_analysis && any_caller_enabled`, gate on `if (!orf_catalogue_active || params.skip_plastid || !params.contrasts)`, delete the `else if` at `:254` (making `:253` a plain `if`), and remove the now-unused `novel_source_configured` or strict-syntax linting will flag it.

7. **Fix the plastid warning's *guard*, not its wording** (`:201`). The wording at `:202` ("the ORF catalogue will still be built") becomes *more* correct after this change. The bug is the `novel_source_configured` conjunct at `:201`, which suppresses the warning for `--extended_orf_analysis true --skip_stringtie true --skip_plastid true` — the exact combination that started this investigation, and one `docs/usage.md:664` promises a warning for. Drop the conjunct.

8. **Update the `:194-200` warning text.** It currently says the flag "has no effect; ORF callers will run against the canonical GTF as usual". After this change caller routing is still a no-op but the catalogue *is* built, so the warning must scope itself to caller annotation and not imply the catalogue is skipped. Do **not** add "and the catalogue will not be built" — that was Revision 1's proposal and it becomes false.

9. **Add a warning when `params.orf_min_callers > enabled_orf_callers.size()`** in `validateInputParameters()`. The schema recommends `2`; a single-caller run then publishes an empty consensus view with no indication.

10. **Update the schema and docs.**
    - `nextflow_schema.json:805-807` — `extended_orf_analysis`'s `description` ("Requires `--skip_stringtie false` or `--novel_gtf <path>`") and `help_text` ("When true but no novel-transcript source is configured, the flag is a no-op") both become false. This is what `nextflow run --help` and the nf-core parameter page render, so it matters more than prose docs.
    - `docs/usage.md:625` — the "no-op" claim is now only true of caller routing.
    - `docs/usage.md:439` — RNA denominator provenance changes in the fallback (relevant if step 5 lands).
    - Add notes on: coordinate-only comparability (Behavior 8), the single-caller consensus caveat, and the cross-caller annotation asymmetry.
    - `docs/output.md:493`, `docs/usage.md:645`, `:654` already describe the new behaviour — **no edit needed**.
    - Stale gating comments: `subworkflows/local/quantify_orf_psite/main.nf:12`, `orf_caller_dispatch/main.nf:43-44` ("empty and unread otherwise"), `main.nf:870`.
    - `CHANGELOG.md` — one `[#NNN]` bullet under `## v1.3.0dev`, per repo convention.

11. **Add the missing test.** No existing test covers extended mode without a novel source: `tests/dotseq.nf.test:15` and `tests/stringtie_extended.nf.test:15` set `skip_stringtie = false`; `tests/novel_gtf.nf.test:15` supplies `novel_gtf`. Add a case with `extended_orf_analysis = true`, `skip_stringtie = true`, no `novel_gtf`. `-profile test` already supplies `contrasts` (`conf/test.config:27`) and is single-caller, so no caller override is needed — and **do not drop `contrasts` to make the test cheaper**, since that is the path that crashed in Revision 1. Follow the `stable_name`/`stable_path` + `tests/.nftignore_orf` pattern from `tests/novel_gtf.nf.test:28-46`. Assert non-zero sizes, not mere existence — every silent failure mode here yields empty or header-only files. Also fix the stale comment at `tests/dotseq.nf.test:20-24` ("Default ribotish + ribocode is enough" — ribocode is off in the test profile).

## Integration

- **Reads**: `params.extended_orf_analysis`, `skip_stringtie`, `novel_gtf`, `skip_plastid`, `contrasts`, `orf_min_callers`, the five caller flags, `ch_gtf`, `ch_hybrid_gtf`.
- **Writes**: new `orf_catalogue/` and `orf_quantification/` in configurations that previously produced neither. Runs that already had a novel source are unaffected — same `ch_full_hybrid_gtf`, same outputs.
- **Order**: predicates must precede `:421`; hence step 1 before step 3.
- **Downstream**: `QUANTIFY_ORF_PSITE` consumes `catalogue_bed12`; ORF-DTE consumes `catalogue_tsv` and the ORF count matrix. Neither needs modification.
- **Snapshots**: all three extended tests configure a novel source, and `default.nf.test` leaves the flag false, so no existing snapshot should move. `BUILD_FULL_HYBRID_GTF` running more often adds no published files (`conf/modules.config:718`). Verify rather than assume. No `tests/tags.yml` change needed — it globs `**.nf.test`.
- **Process names** must not be renamed: they appear in `nf_core_riboseq_software_mqc_versions.yml`, which every pipeline snapshot captures. Tempting for `*_FULL_HYBRID` in fallback runs; don't.

## Assumptions

Fixed rules:

1. `ch_hybrid_gtf` is always defined (`:311`), so `BUILD_FULL_HYBRID_GTF` runs without a novel source.
2. The catalogue annotation must remain a superset of every enabled caller's annotation. Satisfied by construction via the awk — this is why Alternative A was chosen over a `ch_gtf` fallback.
3. Caller routing continues to require a real novel source; there is no hybrid transcriptome without one.
4. `--extended_orf_analysis` remains the catalogue switch. This change does **not** enable the catalogue in default runs.
5. Running `BUILD_FULL_HYBRID_GTF` in fallback runs is acceptable: one gawk streaming pass, nothing published without `--save_reference`.

Configurable behaviour (unchanged): `--skip_plastid` gates ORF quantification and ORF-DTE; `--contrasts` gates DTE; `--skip_orf_collapse` controls smORF collapse; `--orf_min_callers`/`--orf_min_samples` control the consensus view.

## Validation

1. **The changed path works** — `--extended_orf_analysis true`, `--skip_stringtie true`, no `--novel_gtf`, `-profile test` (which supplies contrasts and is single-caller), `--skip_plastid false`. Expect: `BUILD_FULL_HYBRID_GTF` runs, catalogue processes execute, `orf_quantification/` is created, and **no** `EXTENDED_ORF_SECOND_PASS_ALIGN` tasks. Assert `orf_catalogue/*.bed12` and `orf_quantification/*.tsv` are non-empty. This is the configuration that crashed under Revision 1 — it must be the first check, not the fourth.
2. **No regression with a novel source** — `tests/stringtie_extended.nf.test` passes with an unchanged snapshot.
3. **Fallback annotation equals the reference** — with no novel source, assert `BUILD_FULL_HYBRID_GTF`'s output is identical to the input `--gtf` (run with `--save_reference` and diff). Replaces Revision 1's vague "fallback equivalence" prose.
4. **Catalogue comparability on coordinates** — two runs differing only by `--novel_gtf`; compare `orf_catalogue.tsv` on `chrom/start/end/strand` tuples, **not** on `orf_id` (Behavior 8).
5. **`--canonical_gtf` superset holds** — run with a user-supplied `--canonical_gtf` whose transcript IDs are absent from `--gtf`, no novel source. Assert catalogue rows keep non-empty `gene_id` and that multi-exon ORFs retain multi-block BED12. This is the case that made Revision 1's `ch_gtf` fallback unsafe; it must be proven safe under Alternative A.
6. **dotseq no longer over-validates** — `--translational_efficiency_method dotseq --extended_orf_analysis true --skip_stringtie true` with contrasts and plastid. Expect the run to proceed, and specifically **not** to error with an empty requirements list.
7. **Default runs unchanged** — `tests/default.nf.test` passes with an unchanged snapshot and no catalogue.
8. **Warnings fire** — `:194` warning still fires and no longer claims the catalogue is skipped; the plastid warning fires for `--skip_plastid true` (currently suppressed); the new `--orf_min_callers 2` + single-caller warning fires.
9. **No redundant hybrid quant** (if step 5 lands) — assert `quantification/salmon_hybrid/` is absent in the new test.
10. **Linting** — `nf-core lint` (schema text changed) and `nextflow lint` (removed variables, new predicate).

## Questions or ambiguities

- **(Resolved)** Fallback annotation: not `ch_gtf` via a ternary, but `ch_full_hybrid_gtf` built unconditionally — chosen because it satisfies the superset invariant even with a user-supplied `--canonical_gtf`.
- **(Open, non-blocking)** Should the catalogue be available in default runs, whenever ≥1 caller is enabled without `--extended_orf_analysis`? Assumed **no** — it would change default outputs for every user. Both reviewers suggested a dedicated `--orf_catalogue` flag (defaulting to `extended_orf_analysis`) as the clean long-term answer, since that flag now conflates two meanings at the user interface. Out of scope.
- **(Open, non-blocking)** Should the ORF-DTE path gain an RNA-seq presence check? Pre-existing gap, widened here.

## Self-Review

**What changed from Revision 1 and why.** Revision 1 added a `ch_catalogue_gtf = extended_orf_active ? ch_full_hybrid_gtf : ch_gtf` ternary and left `BUILD_FULL_HYBRID_GTF` gated. Both reviewers independently found two fatal problems: `:710`/`:717` read `BUILD_FULL_HYBRID_GTF.out.output` directly, so opening the `:699` block with the process un-invoked crashes at construction; and `:725`/`:751` would receive `channel.empty()`, silently producing zero tasks — the very failure mode Revision 1 claimed to prevent. Revision 1's own step 4 ("audit every use of `ch_full_hybrid_gtf`") would not have found either, because they are uses of the process output object, not the channel.

Building the GTF unconditionally removes both problems, deletes two steps, and closes a third hole: the `ch_gtf` fallback was unsafe with a user-supplied `--canonical_gtf` (the schema recommends MANE Select for human, whose IDs are absent from an Ensembl GTF), where `orfnormalise.py` degrades silently into intron-spanning single-block ORFs and wrong peptides.

**Logic** — the hoist in step 1 is a genuine prerequisite: `enabled_orf_callers` at `:450` currently sits after the `:421` guard that must consume it. Neither review flagged this as blocking; it was found while verifying Alternative A and would have produced an immediate "no such variable" error.

**Corrections to Revision 1's factual claims** — its Context cited `orf_agreement_min_callers` as the adaptive single-caller threshold; that variable is dead code and the live threshold is `params.orf_min_callers`. Its Goal claimed output "directly comparable" to a novel-source run; ORF IDs are positional so comparability is coordinate-only. Its Efficiency section claimed "no runtime cost", which ignored a full STAR genome index build. Its open question about `#204` shifting line numbers was stale — both reviewers verified every line still resolves.

**Reviewer disagreement, resolved.** A said `BUILD_FULL_HYBRID_GTF` publishes only under `--save_reference`; B said it publishes to `quantification/salmon_hybrid`. Both are right: the path is `:716`, the `--save_reference` condition is `:718`. Consequence: no snapshot impact, and the misleading path is visible only with `--save_reference`. A and B also differed on whether `ch_gtf`'s channel type needed a defensive `.first()`; Alternative A dissolves the question, since no new channel is introduced.

**Edge cases** — added the Ribo-seq-only, `--orf_min_callers 2`, `--canonical_gtf` and cross-caller-asymmetry cases, none of which Revision 1 covered.

**Remaining risks** — (a) the predicate is still duplicated between `main.nf` and `dotseqPrerequisitesError`; a shared helper is impractical because validators run before the workflow and cannot import a workflow-scope `def` without a wider refactor, so this is accepted and noted. (b) Step 5 is the largest single piece. It is confirmed in scope, so the PR carries both a gating change and a denominator change; if review pushes back on size, it is the separable half. (c) Validation 5 (`--canonical_gtf`) requires constructing a test annotation with IDs absent from the reference, which is fiddly; if it is skipped, the superset invariant rests on reading the awk rather than on a test.

## Revision History

| Rev | Date | Change |
|---|---|---|
| 1 | 2026-07-31 | Initial plan: second predicate + `ch_catalogue_gtf` ternary with a `ch_gtf` fallback. |
| 2 | 2026-07-31 | Rewritten after dual review. Adopts Alternative A (build the full hybrid GTF unconditionally), dropping the ternary. Adds the mandatory `enabled_orf_callers` hoist, fixes the dotseq gate and plastid warning guard as *conditions* rather than messages, adds schema updates, corrects the dead-code and comparability claims, and adds five edge cases and four validations. |
