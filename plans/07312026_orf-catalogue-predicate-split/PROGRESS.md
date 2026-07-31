# Progress: Decouple the ORF catalogue from the novel-transcript source

**Plan:** [PLAN.md](PLAN.md)
**Created:** 2026-07-31
**Status:** Implemented in #222

## Pipeline

| Step | Status | Notes |
|---|---|---|
| Plan written | ✅ Done | Rev 1, 2026-07-31 |
| Agent review (`plan-reviewer`) | ✅ Done | Dual: `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` |
| Plan revised | ✅ Done | Rev 2 — adopts Alternative A; 3 Criticals resolved |
| Manual review | ✅ Done | Rev 2 approved; implemented in #222 |
| Implementation | ✅ Done | #222, branch `feat/orf-catalogue-without-novel-source` |
| Verification | ⬜ Not started | `code-reviewer` + `plan-manager` after implementation |

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
