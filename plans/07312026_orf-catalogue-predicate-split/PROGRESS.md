# Progress: Decouple the ORF catalogue from the novel-transcript source

**Plan:** [PLAN.md](PLAN.md)
**Created:** 2026-07-31
**Status:** Implemented and reviewed in #222 — snapshot outstanding

## Pipeline

| Step | Status | Notes |
|---|---|---|
| Plan written | ✅ Done | Rev 1, 2026-07-31 |
| Agent review (`plan-reviewer`) | ✅ Done | Dual: `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` |
| Plan revised | ✅ Done | Rev 2 — adopts Alternative A; 3 Criticals resolved |
| Manual review | ✅ Done | Rev 2 approved; implemented in #222 |
| Implementation | ✅ Done | #222, branch `feat/orf-catalogue-without-novel-source` |
| Verification | ✅ Done | Dual `CODE_REVIEW_A/B.md` + `COVERAGE.md`; fixes in `4f06519` |

## Context

Arose from a Ribo-seq run where the ORF catalogue and `orf_quantification/` were
absent. Diagnosis took three passes and two corrections:

1. First read: catalogue skip was silent, no warning existed. **Wrong** — the
   warning exists at `utils_nfcore_riboseq_pipeline/main.nf:199`; the initial
   grep covered only `workflows/riboseq/main.nf`.
2. Second read: fail-fast on the contradictory flag combination. **Wrong** —
   `docs/usage.md:625` documents the no-op as deliberate, so users can compose
   flags incrementally.
3. Third read: `--skip_stringtie true` left `extended_orf_active` false, which
   skipped the catalogue (`:477`); `--skip_plastid true` independently removed
   `orf_quantification/` (`:546`). Confirmed when a rerun with both set to
   `false` started `EXTENDED_ORF_SECOND_PASS_ALIGN`.

The behaviour was working as designed and as documented. The plan is therefore a
**feature change**, not a bug fix: the user wants the catalogue to follow
`--extended_orf_analysis` alone, independent of StringTie.

## Decisions

| Date | Decision |
|---|---|
| 2026-07-31 | ~~Fallback annotation is `ch_gtf`, not `ch_canonical_gtf`~~ — **superseded by Rev 2**: no fallback channel at all |
| 2026-07-31 | Caller routing keeps requiring a real novel source; only the catalogue chain decouples |
| 2026-07-31 | Add a second predicate rather than redefining `extended_orf_active` |
| 2026-07-31 | Dropped the earlier warning/docs wording tweaks — made moot by this change |
| 2026-07-31 | **Rev 2:** build `BUILD_FULL_HYBRID_GTF` unconditionally instead of adding a `ch_gtf` fallback ternary — Rev 1 crashed at workflow construction |
| 2026-07-31 | Hoist `enabled_orf_callers` above `:421` — prerequisite, not cleanup |
| 2026-07-31 | Fix the dotseq gate and plastid warning as *conditions*, not messages |

## Review outcome (Rev 1)

Both reviewers independently concluded **do not implement as written**. Three Criticals,
all now resolved in Rev 2:

1. `main.nf:710`/`:717` read `BUILD_FULL_HYBRID_GTF.out.output` directly, so Rev 1's
   gating crashed at construction. A grep for `ch_full_hybrid_gtf` — which Rev 1's own
   step 4 prescribed — does not find those lines.
2. The dotseq fix targeted the message, not the gate; as written it would error with an
   empty requirements list.
3. The `ch_gtf` fallback was unsafe with a user-supplied `--canonical_gtf`, silently
   producing intron-spanning ORFs and wrong peptides.

Reviewer disagreement on `BUILD_FULL_HYBRID_GTF`'s publishDir resolved: both partly
right — path is `conf/modules.config:716`, `--save_reference` condition is `:718`. No
snapshot impact.

## Open items

- Should the catalogue also run in default (non-extended) runs? Assumed no.
- ~~`#204` may have shifted the `:699` region~~ — resolved: both reviewers verified every line number still resolves at `f4a3753`.
- No existing test covers extended mode without a novel source — the plan adds one.

## Verification outcome

Coverage verdict **INCOMPLETE — 9 gaps**: steps 1-9 (all code) DONE with no deviations;
every gap in step 10/11 sub-items or unrun validations. Four doc/comment misses fixed in
`4f06519`.

Both code reviewers independently found:

1. The new test's `!salmon_hybrid.exists()` assertion was **false by construction** — Nextflow
   creates a publishDir target directory whenever the process runs, even when `saveAs` returns
   null. Confirmed from committed snapshots (`default.nf.test.snap` genome dirs,
   `dotseq.nf.test.snap:1621` salmon_hybrid/star_index, both with no files beneath). The test
   would have gone red in CI. Fixed to assert on contents.
2. **The RNA-denominator reuse reintroduces the gene-id namespace mismatch** that building the
   hybrid GTF unconditionally was chosen to avoid — one stage downstream. Catalogue genes resolve
   against `ch_gtf ∪ canonical-only`; the reused primary matrix derives from `ch_gtf` alone, so a
   user-supplied `--canonical_gtf` puts canonical-only genes' ORFs outside the denominator and
   `DTE_COUNTS_PREP` drops them with only a stderr count. Documented; a code-level signal is still
   an open decision.

Reviewer disagreement, resolved: A said the publishDir directory *is* created (test fails), B said
it is not (test passes for the wrong reason). Checked the committed snapshots — **A is correct**.
Both agreed on the fix regardless.

Corrections to earlier artifacts: the plan's `as boolean` justification was wrong (Groovy `&&`
already yields a Boolean); the plan's "no snapshot impact" claim was wrong (an empty directory does
enter `stable_name`); and the plan-review premise that the dotseq validator "would error with an
empty requirements list" was wrong — the old `if/else if` pair covered both branches. The resulting
code is correct in all three cases.

## Outstanding

- **Blocking:** record and commit `tests/orf_catalogue_no_novel_source.nf.test.snap`. Needs a real
  run; not possible locally (arm64 emulation, seqkit crashes under it).
- Validations 2, 3, 4, 5, 7, 10 unrun. Validation 5 (`--canonical_gtf` superset) is the one the
  plan flagged as load-bearing.
- Open decision: how to signal the `--canonical_gtf` denominator mismatch — docs only (done),
  a `log.warn`, or gate the matrix reuse on `!params.canonical_gtf`.
- Deferred: single `enabledOrfCallers()` helper to collapse three encodings of the caller set.
