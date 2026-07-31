# Plan review — Decouple the ORF catalogue from the novel-transcript source (Reviewer B)

**Plan**: `plans/07312026_orf-catalogue-predicate-split/PLAN.md`
**Repo state reviewed**: `HEAD = 6b17ac7` (= `origin/dev` `f4a3753` + a CHANGELOG-only commit). All `workflows/`, `subworkflows/`, `docs/` line numbers quoted in the plan are therefore **current** — the plan's open question about `#204` having shifted the `:699` region is resolved: nothing moved.

**Verdict**: the goal is right and the diagnosis of the two-purposes-one-predicate problem is accurate. But the implementation outline, if followed literally, **fails at workflow construction** in the exact configuration it targets, and two of its factual premises are wrong. Do not implement as written.

---

## 1. Verification of the plan's factual claims

| Plan claim | Status |
|---|---|
| `:361-362` predicate pair, `:366`/`:421`/`:439` caller routing, `:477`/`:546`/`:699` catalogue gates | **Correct**, all verified verbatim |
| `enabled_orf_callers` built at `:450-456`; `orf_agreement_min_callers` at `:463` | Lines correct, **but see 4.2 — that variable is dead code** |
| `BUILD_FULL_HYBRID_GTF` = full reference + novel genes appended | **Correct and stronger than stated**: `assets/build_full_hybrid_gtf.awk:13-26` prints file 1 (`ch_gtf`) verbatim, then only rows from file 2 whose `transcript_id` (or `gene_id`) is absent. With zero novel entries the output is byte-identical to `ch_gtf` |
| `dotseqPrerequisitesError()` at `utils_nfcore_riboseq_pipeline/main.nf:244-260` duplicates the predicate | **Correct** (`:247-248`) |
| `extended_orf_active` used at seven sites in `main.nf` | **Correct** (`:362` def + `:366`, `:421`, `:439`, `:477`, `:546`, `:699`) |
| "again in two validators" | **Loose**: `dotseqPrerequisitesError` recomputes the composite (`:248`); `validateInputParameters` recomputes only `novel_source_configured` (`:197`). Both still need edits |
| No existing test covers extended-ORF without a novel source | **Correct**: `dotseq.nf.test:15` and `stringtie_extended.nf.test:15` set `skip_stringtie=false`; `novel_gtf.nf.test:15` supplies `novel_gtf` |
| `docs/usage.md:625` says the flag "is a no-op" | **Correct**, and incomplete — see 4.3 |
| Existing snapshots unaffected | **Correct.** `default.nf.test` leaves `extended_orf_analysis` false (`nextflow.config:150`); the three extended tests keep a novel source. No test asserts on `workflow.stdout`, so the warning-text rewrites in steps 6-7 are snapshot-safe |
| `ch_gtf` is the correct fallback because full_hybrid == ch_gtf when there are no novel genes | **True only when `--canonical_gtf` is unset** — see 3.2 |
| "no runtime cost" (Self-Review, Efficiency) | **Wrong** — see 5 |

---

## 2. Critical logic flaw: the plan crashes the pipeline it is trying to enable

`workflows/riboseq/main.nf:709-718` references the **process output object directly**, not the `ch_full_hybrid_gtf` channel:

```groovy
GFFREAD_FULL_HYBRID(
    BUILD_FULL_HYBRID_GTF.out.output,          // :710
    ch_fasta
)
...
STAR_GENOMEGENERATE_FULL_HYBRID(
    ch_fasta.map { ... },
    BUILD_FULL_HYBRID_GTF.out.output.map { ... }  // :717
)
```

Plan step 2 says *"Leave `BUILD_FULL_HYBRID_GTF` and its `if (extended_orf_active)` guard alone"*. Plan step 3 changes `:699` to `if (orf_catalogue_active && !params.skip_plastid)`. Those two decisions are incompatible: with `--extended_orf_analysis true --skip_stringtie true`, `extended_orf_active` is false so the process at `:421-428` is never invoked, but the block at `:699` is now entered, and Nextflow raises at construction time:

```
Access to 'BUILD_FULL_HYBRID_GTF.out' is undefined since the process
'BUILD_FULL_HYBRID_GTF' has not been invoked before accessing the output attribute
```

I confirmed this is a hard error, not a warning, with a minimal reproducer on Nextflow 25.10.4 (process guarded by a false `if`, `.out` read inside a true `if`) — it fails before any task is submitted.

This is not a corner case. `conf/test.config:27` sets `contrasts` for **every** `-profile test` run, and `skip_plastid` defaults to false (`nextflow.config:125`). So the new nf-test the plan asks for in step 9 — `extended_orf_analysis=true, skip_stringtie=true`, `-profile test` — enters the `if (ch_contrasts_file)` block at `:655`, enters the `:699` block, and dies immediately. Plan Validation 1 would *not* catch it (it doesn't mention contrasts); only Validation 4 would, and by then the implementer has a broken branch and no diagnosis.

Step 4 gestures at this (*"Check whether any other consumer of `ch_full_hybrid_gtf` should also fall back — audit every use"*) but the two dangerous references are **not** uses of `ch_full_hybrid_gtf` — they are uses of `BUILD_FULL_HYBRID_GTF.out`, which a grep for the channel name will not surface. The plan must name them explicitly.

Note also `:725` and `:751` inside the same block do read `ch_full_hybrid_gtf`, which under the plan's design would be `channel.empty()` — so even if the `.out` references were fixed, `FASTQ_ALIGN_STAR_FULL_HYBRID` and `QUANTIFY_HYBRID_RNA` would be constructed with an empty GTF channel and silently produce zero tasks. That is precisely the failure mode the plan's own Self-Review says step 2 exists to prevent; the plan closed it at `:487` and left it open at `:725`/`:751`.

---

## 3. Assumptions

### 3.1 Correct and worth keeping

- Assumption 2 (caller routing still needs a real novel source) is right: `orf_caller_dispatch/main.nf:210-228` swaps in `ch_hybrid_transcriptome_bam`, which only exists under `if (extended_orf_active)` at `main.nf:366-394`. Redefining `extended_orf_active` would indeed have enabled a second STAR pass with no hybrid transcriptome.
- Assumption 3 (the catalogue stays behind `--extended_orf_analysis`) is the right call for a default-output-stability reason the plan states correctly.
- The channel-typing question is a non-issue: `orftable_fasta_gtf_buildorfcatalogue/main.nf` applies `.first()` to both the GTF and FASTA inputs internally, so passing `ch_gtf` (a single-item queue in some `PREPARE_GENOME` branches, e.g. `prepare_genome/main.nf:78`) instead of the `.first()`-ed `ch_full_hybrid_gtf` is safe. Worth stating in the plan so nobody "fixes" it later.

### 3.2 Assumption 1 is false when `--canonical_gtf` is supplied

The plan asserts `ch_gtf` and full-hybrid-with-no-novel are "equivalent". They are equivalent only because `ch_hybrid_gtf` defaults to `ch_canonical_gtf` (`main.nf:311`) and `ch_canonical_gtf` is normally derived *from* `ch_gtf` by AGAT longest-isoform extraction (`prepare_genome/main.nf:141-149`), so its transcript IDs are a subset.

But `--canonical_gtf` lets the user supply an arbitrary file (`prepare_genome/main.nf:131-136`), and the schema **recommends MANE Select for human** (`nextflow_schema.json`, `canonical_gtf.help_text`). MANE RefSeq IDs (`NM_…`) do not exist in an Ensembl `--gtf`. In that configuration:

- **Today, extended mode**: full_hybrid = `ch_gtf` + every hybrid row whose transcript is absent from `ch_gtf` — so the MANE transcripts **are** carried into the catalogue annotation (`build_full_hybrid_gtf.awk:24`).
- **Under the plan, fallback mode**: the catalogue gets plain `ch_gtf`, which contains **none** of the transcripts Ribo-TISH and Ribotricer just called ORFs against (`orf_caller_dispatch/main.nf:54-55`, `:133-135`).

The catalogue normaliser degrades silently rather than failing (`modules/nf-core/custom/orfnormalise/templates/orfnormalise.py`):

- ribotish, `:557-568`: `tx is None` → the ORF becomes a **single-block BED12** over the whole genomic span, so any spliced ORF gets intronic sequence baked into the catalogue FASTA and the P-site BED.
- ribotricer, `:599-616`: same single-span fallback.
- ribocode, `:473-491`: no transcript and no `ORF_gstart` → `continue`, the ORF is **dropped**.
- rpbp, `:715` / price, `:806-819`: `gene_id` becomes `""` or a concatenated multi-gene string → the `orf_to_gene` join in `DTE_COUNTS_PREP` drops those rows.

So the plan's fallback can produce a wrong-but-plausible catalogue for `--canonical_gtf` users, with no error. This is the single strongest argument for alternative A below.

### 3.3 Unstated assumption: at least one RNA-seq sample exists

`main.nf:316-322` only enforces "RNA-seq BAMs required" when `run_stringtie` is true. The plan's new path is `skip_stringtie=true`, so nothing checks it, and `ch_rnaseq_reads` (`:720`) can be empty → `DTE_COUNTS_PREP` never fires → no ORF-DTE output, no message. Pre-existing (reachable today via `--novel_gtf`), but the plan widens the door. Worth one sentence in the plan and, ideally, a check.

---

## 4. Gaps — things the plan does not mention that it must

### 4.1 Steps 5 and 7 are under-specified in a way that leaves the bug in place

**Step 5 (`dotseqPrerequisitesError`)**: the plan says the *"a novel-transcript source"* clause "must go". Deleting only that clause (`utils_nfcore_riboseq_pipeline/main.nf:254`) leaves the guard at `:251` reading `!extended_orf_active`, which is still true for the newly valid combination — the pipeline errors with an **empty** requirements list: `"...requires: . Pick anota2seq or deltate..."`. The condition at `:248` must change to drop `novel_source_configured`, not just the message line.

**Step 7 (plastid warning)**: the plan says update it *"if its wording implies…"*. The wording (`:202`, "the ORF catalogue will still be built") is already correct post-change. The **condition** is the bug: `:201` requires `novel_source_configured`, so the user who runs `--extended_orf_analysis true --skip_stringtie true --skip_plastid true` — the combination the plan's Behavior/edge-cases section calls out as "the combination that produced the original confusion" — gets **no warning at all**. Drop `novel_source_configured` from `:201`.

### 4.2 The single-caller edge case rests on dead code, and the real trap is different

Plan Context bullet 2 and the "Single caller" edge case both rely on `orf_agreement_min_callers` (`main.nf:463`) adapting to the caller count. That variable — along with `rank_aggregation_callers` (`:460`), `ch_enabled_orf_callers` (`:466`) and `ch_rank_aggregation_callers` (`:467`) — is **defined and never consumed anywhere in the pipeline**. The threshold the merger actually uses is the static user parameter: `conf/modules.config:1081` passes `--min-callers ${params.orf_min_callers}`, default `1` (`nextflow.config:164`).

Consequences:

- The plan's reasoning ("threshold becomes 1, every ORF trivially meets it") reaches the right answer for the wrong reason. Fine at defaults, but it means nothing adapts.
- The real trap: `nextflow_schema.json` (`orf_min_callers.help_text`) actively recommends `--orf_min_callers 2` "for a higher-confidence catalogue". Combine that with a single-caller run — which this change makes *much* easier to reach, and which is exactly what `-profile test` produces (`conf/test.config:33` sets `skip_ribocode=true`, leaving Ribo-TISH alone) — and `orf_catalogue/consensus/*` is published **empty apart from headers**, silently. That deserves a warning when `params.orf_min_callers > enabled_orf_callers.size()`, not a docs footnote.

### 4.3 Documentation and housekeeping the plan omits

- **`nextflow_schema.json:805-806`** — the `extended_orf_analysis` `description` says *"Requires --skip_stringtie false or --novel_gtf <path>"* and the `help_text` says *"When true but no novel-transcript source is configured, the flag is a no-op"*. Both become false. This matters more than `docs/usage.md` because it is what `nextflow run --help` and the nf-core website parameter page render. **Not mentioned in the plan at all.**
- **`docs/usage.md:439`** — "The RNA-seq matrix used here is quantified against the full reference transcriptome augmented with the novel intergenic transcripts" is wrong in the fallback whatever design is chosen (see 5).
- **`docs/usage.md:431`** — pre-existing inaccuracy: it claims ORF-level DTE requires `--te_quantification_method plastid_psite`, but the code gate (`main.nf:699`) never reads that param. Cheap to fix while in the area.
- **`docs/output.md:493`** already reads *"Produced only when `--extended_orf_analysis true` is set and at least one ORF caller is enabled"* — i.e. the docs already describe the post-change behaviour. Good supporting evidence that the coupling is a bug rather than a design choice; no edit needed there.
- **`CHANGELOG.md`** — the repo keeps a per-PR entry under `## v1.3.0dev` (see `:10-25`). Not mentioned.
- **`subworkflows/local/quantify_orf_psite/main.nf:12`** and **`orf_caller_dispatch/main.nf:43-44`** ("empty and unread otherwise") carry gating comments that go stale under either design.

---

## 5. Efficiency — the Self-Review's "no runtime cost" is wrong

Gating workflow *construction* is free, yes. But the plan newly enters the `:699` block, and that block is expensive and, without a novel source, **entirely redundant**:

- `GFFREAD_FULL_HYBRID` (`:709`) rebuilds a transcript FASTA from a GTF that equals `ch_gtf`.
- `STAR_GENOMEGENERATE_FULL_HYBRID` (`:715`) builds a second full STAR index — tens of GB of RAM and disk on a mammalian genome.
- `FASTQ_ALIGN_STAR_FULL_HYBRID` (`:722`) re-aligns **every RNA-seq library** against it.
- `QUANTIFY_HYBRID_RNA` (`:746`) re-quantifies them with Salmon.

The comment at `:700-708` states the rationale precisely: this exists *because* novel genes have no row in the primary Salmon matrix. With no novel genes that rationale evaporates — the pipeline would pay for a second index and a second alignment pass to recompute a matrix it already has at `:588` / `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled`.

And the substitution is already supported: `modules/local/dte_counts_prep/templates/dte_counts_prep.py:17-21` and `:126-137` explicitly handle *"a Salmon all-sample matrix carrying its Ribo-seq columns alongside the wanted RNA-seq ones"* by dropping the overlapping primary-role columns. Both matrices are `counts_gene_length_scaled` from the same subworkflow, so the shapes are identical and the gene IDs both derive from `ch_gtf`. Reusing the primary matrix in the fallback is a drop-in that is *also* scientifically preferable: the ORF-level RNA denominator becomes literally the same numbers as the gene-level one.

---

## 6. Scientific correctness and downstream comparability

### 6.1 Cross-caller annotation scope in the fallback (biases consensus, unchanged from today but newly prominent)

In the fallback, callers do **not** share an annotation (`orf_caller_dispatch/main.nf`):

| Caller | Fallback annotation | Extended annotation |
|---|---|---|
| Ribo-TISH `predict` | `ch_canonical_gtf`, no `-a` (`:55`, `:96`) | hybrid `-g` + canonical `-a` (`:69-73`) |
| Ribotricer | `ch_canonical_gtf` (`:135`) | `ch_hybrid_gtf` (`:134`) |
| Rp-Bp / PRICE | `ch_gtf` multi-isoform (`:65`) | `ch_full_hybrid_gtf` (`:64`) |
| RiboCode | `ch_gtf` (`:223`) | `ch_hybrid_gtf` (`:222`) |

So one-transcript-per-gene callers are pooled with multi-isoform callers. An isoform-specific ORF is invisible to Ribo-TISH/Ribotricer by construction, so it can never reach 2-of-2 agreement in the default Ribo-TISH + RiboCode pair — `called_by_*` and any `--orf_min_callers 2` consensus systematically penalise exactly the ORFs the multi-isoform annotation exists to find. This asymmetry is **pre-existing** (extended mode has the same split) and the plan is right not to touch caller routing. But the plan's headline is "the catalogue is now generally available", so this becomes the first thing a careful user hits, and the plan's docs step mentions only the single-caller case. It belongs in the docs note.

One thing the plan gets right by omission: Ribo-TISH ORF *classification* is comparable across modes, because `-a` is the canonical backbone in extended mode and `-g` is the canonical backbone in fallback mode — the class labels come off the same annotation either way.

### 6.2 The Goal's comparability claim is too strong

Plan line 9: *"its output is directly comparable to a novel-source run minus the novel genes."* Catalogue ORF IDs are **positional**: `modules/nf-core/custom/orfmerge/templates/orfmerge.py:206-234` sorts clusters by coordinate and assigns `orf_{idx+1:08d}`. Inserting novel-gene ORFs renumbers every ORF downstream of them. Two runs are comparable by **coordinates**, never by ORF ID — and the ORF-level DTE tables under `dte/orf_level/` are keyed by those IDs, so a user joining two runs' DTE results on `orf_id` gets a silent, systematic mis-join. Add to that the collapse step: `MMSEQS_EASYCLUSTER` + `CUSTOM_ORFCOLLAPSE` cluster catalogue peptides cohort-globally, so adding novel peptides can change cluster membership and representative selection for ORFs the two runs share.

The plan should reword the Goal to "comparable on coordinates" and say so in the docs. Validation 3 must diff on `chrom/start/end/strand`, not on `orf_id`.

---

## 7. Validation sufficiency

| Plan step | Assessment |
|---|---|
| 1. Changed path works | **Under-specified.** Add `--contrasts` explicitly (or state that `-profile test` supplies it) — without contrasts it passes while the Critical bug is live. Also assert `orf_catalogue/orf_catalogue.bed12` and `orf_quantification/orf_psite_counts.tsv` are **non-empty**, not just present: every silent-failure mode in this change produces empty-or-header-only files |
| 2. No regression with a novel source | Adequate |
| 3. Fallback equivalence | Weakest step, as the plan admits. Make it concrete: two runs identical except `--novel_gtf`, then compare `orf_catalogue.tsv` **on coordinate tuples** (see 6.2). Also assert `BUILD_FULL_HYBRID_GTF` did not run (or, under alternative A, that its output is identical to `--gtf`) |
| 4. dotseq no longer over-validates | Good, and it is the only step that catches the Critical bug — reorder it earlier, or fold contrasts into step 1 |
| 5. Default runs unchanged | Adequate |
| 6. Warning accuracy | Extend: assert the **plastid** warning fires for `extended_orf_analysis=true, skip_stringtie=true, skip_plastid=true` (currently suppressed, 4.1) |
| — | **Missing**: a case with `--canonical_gtf` set and no novel source, asserting catalogue rows keep non-empty `gene_id` and multi-block BED12 (3.2) |
| — | **Missing**: `--orf_min_callers 2` with one caller → assert the consensus view is empty *and* that a warning was emitted (4.2) |

CI note: the new nf-test in step 9 is a full pipeline run that now includes a STAR index build and an extra alignment pass. Under alternative B below that cost disappears; under the plan as written it does not. Either way, do **not** drop `contrasts` from the new test to make it cheaper — that is the only part of the path that currently crashes. `tests/tags.yml` needs no change (it globs `**.nf.test`), and the test should follow the existing `stable_name` / `stable_path` + `tests/.nftignore_orf` pattern from `tests/novel_gtf.nf.test:28-46`.

---

## 8. Alternatives

**A. Build the full GTF unconditionally (recommended).** `ch_hybrid_gtf` is *always* defined — it defaults to `ch_canonical_gtf` at `main.nf:311` and is never empty. So `BUILD_FULL_HYBRID_GTF(ch_gtf.combine(ch_hybrid_gtf) …)` works with no novel source: it yields `ch_gtf` plus any canonical-only transcripts. Change `:421` from `if (extended_orf_active)` to `if (extended_orf_active || orf_catalogue_active)` and stop there. This:

- needs **no** `ch_catalogue_gtf` ternary and no change at `:487` — one channel, one truth;
- fixes the Critical bug at `:710`/`:717` and the empty-channel problem at `:725`/`:751` for free;
- **preserves superset-ness even with `--canonical_gtf`** (3.2), because the canonical-only transcripts get appended exactly as they do today in extended mode;
- costs one `gawk` task.

Cost/caveats to handle: the module publishes to `${params.outdir}/quantification/salmon_hybrid` (`conf/modules.config:712-716`), so a fallback run would emit a `full_hybrid_reference.gtf` under a "hybrid"/"salmon" path with nothing hybrid in it — either disable that `publishDir` when `!extended_orf_active`, or accept it and document it. Do **not** rename the process alias to something neutral, tempting as it is: process names appear in `nf_core_riboseq_software_mqc_versions.yml`, which every pipeline snapshot captures.

**B. Keep the plan's ternary, and skip the hybrid-RNA machinery in the fallback.** If A is rejected for output-naming reasons, then step 3's edit to `:699` must be paired with an explicit fallback inside the block: when `!extended_orf_active`, skip `GFFREAD_FULL_HYBRID`, `STAR_GENOMEGENERATE_FULL_HYBRID`, `FASTQ_ALIGN_STAR_FULL_HYBRID` and `QUANTIFY_HYBRID_RNA` entirely and feed `DTE_COUNTS_PREP` (`:762`) with `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled` (already column-compatible, see 5). This is the cheapest correct resolution of the Critical bug and removes a redundant STAR index build. It does **not** fix 3.2 — A does.

Best combination: **A for the catalogue annotation + B for the RNA denominator.** They are independent and both are net simplifications.

**C. Cosmetic but worth it: define the predicates together.** The plan puts `orf_catalogue_active` at `:456`, ~95 lines from `extended_orf_active` at `:362`. `enabled_orf_callers` (`:450-456`) is pure param logic with no channel dependency — hoist it to just above `:361` and define all three predicates in one block. A maintainer scanning for "what gates what" then reads one place instead of three. (The plan's Remaining-risks note about deferring a shared helper for the validators is the right call — the validators run before the workflow and cannot import a workflow-scope def without a refactor that widens the diff for no behavioural gain.)

---

## 9. Action items

### Critical

1. **Fix the `BUILD_FULL_HYBRID_GTF.out` references.** `main.nf:710` and `:717` read the process output object directly; with the plan's gating the process is never invoked and Nextflow fails at construction (verified empirically). Adopt alternative A (build the GTF unconditionally) or alternative B (skip the hybrid-RNA block and use the primary Salmon matrix). Note that `:725` and `:751` would also receive `channel.empty()` under the plan as written.
2. **Fix the `dotseqPrerequisitesError` *condition*, not just its message.** Change `utils_nfcore_riboseq_pipeline/main.nf:248`/`:251` to stop requiring `novel_source_configured`; deleting only the `missing <<` clause at `:254` still errors, with an empty requirements list.
3. **Correct Assumption 1.** `ch_gtf` is not equivalent to full-hybrid-with-no-novel when `--canonical_gtf` is supplied (schema recommends MANE for human). The catalogue normaliser then silently emits single-block BED12s, empty `gene_id`s, or drops ORFs (`orfnormalise.py:473-491`, `:557-568`, `:599-616`, `:715`, `:806-819`). Alternative A removes this class of bug.
4. **Add `--contrasts` to Validation 1** (or reorder Validation 4 first) so the crash path is exercised on the first check rather than the fourth.

### Important

5. **Fix the plastid warning's condition** at `utils_nfcore_riboseq_pipeline/main.nf:201` — drop `novel_source_configured`, otherwise the `skip_plastid` combination the plan calls out gets no warning at all.
6. **Update `nextflow_schema.json:805-806`** (`extended_orf_analysis` `description` + `help_text`). Both state the novel-source requirement and the "no-op" behaviour; both become false. Not in the plan.
7. **Correct the Context claim about `orf_agreement_min_callers`.** It, `rank_aggregation_callers`, `ch_enabled_orf_callers` and `ch_rank_aggregation_callers` (`main.nf:460-467`) are dead code; the merger's threshold is `params.orf_min_callers` (`conf/modules.config:1081`, default 1). Add a warning when `params.orf_min_callers > enabled_orf_callers.size()` — the schema recommends `2`, and a single-caller run then publishes an empty consensus view silently.
8. **Correct the Self-Review "no runtime cost" claim** and the Goal's "directly comparable" claim (comparable on coordinates only — ORF IDs are positional, `orfmerge.py:206-234`; the cohort-global peptide collapse can also reshuffle shared ORFs). Validation 3 must diff on coordinates.
9. **Extend the docs step**: `docs/usage.md:439` (RNA denominator provenance in the fallback), the cross-caller annotation-scope asymmetry from 6.1 (one-transcript-per-gene callers pooled with multi-isoform callers biases `called_by_*`), and the ORF-ID-instability caveat. `docs/output.md:493` already matches the new behaviour and needs no edit.
10. **Add a `CHANGELOG.md` entry** under `## v1.3.0dev` per repo convention.

### Optional

11. Hoist `enabled_orf_callers` and define all three predicates in one block near `main.nf:361` (alternative C).
12. Add an RNA-seq-sample presence check for the ORF-DTE path (3.3) — pre-existing gap, widened by this change.
13. Refresh the now-stale gating comments at `subworkflows/local/quantify_orf_psite/main.nf:12` and `orf_caller_dispatch/main.nf:43-44`.
14. Fix `docs/usage.md:431`'s incorrect claim that ORF-level DTE requires `--te_quantification_method plastid_psite` (the code gate never reads it) while in the area.
15. Resolve the plan's open question about single-caller warnings in favour of a warning rather than docs-only — item 7 makes it necessary anyway.
