# Plan Review A — Decouple the ORF catalogue from the novel-transcript source

Reviewed: `plans/07312026_orf-catalogue-predicate-split/PLAN.md`
Repo state: `HEAD = 6b17ac7` (on top of `f4a3753`). **All line numbers cited in the plan still match this HEAD exactly** — `:361`, `:362`, `:366`, `:420-428`, `:439`, `:450-456`, `:463`, `:477`, `:487`, `:546`, `:699` all resolve to what the plan says they do. The plan's open item about `#204` shifting line numbers is stale: nothing moved. Verified by reading `workflows/riboseq/main.nf`, `subworkflows/local/orf_caller_dispatch/main.nf`, `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf`, `subworkflows/local/prepare_genome/main.nf`, `subworkflows/nf-core/orftable_fasta_gtf_buildorfcatalogue/main.nf`, `assets/build_full_hybrid_gtf.awk`, `modules/nf-core/custom/orfnormalise/templates/orfnormalise.py`, `modules/local/dte_counts_prep/`, `conf/modules.config`, `conf/test.config`, `nextflow_schema.json`, `docs/{usage,output}.md` and all four relevant `tests/*.nf.test`.

Verdict: the diagnosis is right and the two-predicate design is the right shape, but **step 3 as written breaks the pipeline in exactly the configuration the plan exists to enable**, and **step 5 as written can make the dotseq validator worse rather than better**. Two Critical items, six Important, plus a materially simpler alternative that removes the first Critical for free.

---

## 1. Logic review

### Confirmed correct

- **The coupling is real and doubly enforced** (plan Context, line 27). Confirmed: `ch_full_hybrid_gtf = channel.empty()` at `main.nf:420`, and `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE` re-applies `.first()` internally (`subworkflows/nf-core/orftable_fasta_gtf_buildorfcatalogue/main.nf:44`: `CUSTOM_ORFNORMALISE ( ch_normalise_in, ch_gtf.first() )`). `channel.empty().first()` emits nothing, so the process would be constructed and produce **zero tasks with no error**. The plan's Self-Review is right that changing only the `if` conditions reproduces the original symptom invisibly. Good catch; steps 3 and 4 must indeed land together.
- **Not editing `extended_orf_active` in place.** Correct and important. `extended_orf_active` is passed as a `val` into `ORF_CALLER_DISPATCH` (`main.nf:439`) and drives six independent branches there (`orf_caller_dispatch/main.nf:63, 94, 133, 210, 221, 226`). Widening it would have switched RiboCode onto `ch_hybrid_transcriptome_bam` — a `channel.empty()` when the second pass didn't run (`main.nf:364`) — silently zeroing RiboCode. The plan avoids this.
- **`ch_gtf` is the right fallback**, and more strongly than the plan argues. `orf_caller_dispatch/main.nf:63-65` **already implements the identical ternary** for Rp-Bp and PRICE: `extended_orf_active ? full_hybrid : ch_gtf`. So `ch_gtf` is the established non-extended counterpart of `full_hybrid` in this codebase, not a new invention. `assets/build_full_hybrid_gtf.awk` confirms the set relation exactly: output = *every row of `ch_gtf`*, then hybrid rows whose `transcript_id` is absent from `ch_gtf`. With no novel transcripts the two are byte-identical. The `docs/usage.md:387` citation against `ch_canonical_gtf` is accurate.
- **`:546` needs only `orf_catalogue_active`** — it is nested inside `if (!params.skip_plastid)` at `:523`. Correct.
- **Snapshot claim holds.** `tests/dotseq.nf.test:15` and `tests/stringtie_extended.nf.test:15` set `skip_stringtie = false`; `tests/novel_gtf.nf.test:15-16` supplies `novel_gtf` with `skip_stringtie = true`. All three configure a novel source. `tests/default.nf.test` never sets `extended_orf_analysis`. No existing snapshot should move.
- **The docs are already ahead of the code**, which corroborates the plan's premise that the coupling was unintended: `docs/output.md:493` ("Produced only when `--extended_orf_analysis true` is set and at least one ORF caller is enabled"), `docs/usage.md:645`, `docs/usage.md:654` ("gates on `--extended_orf_analysis true` plus a non-empty enabled-caller set") and `docs/usage.md:431` all describe the *post-change* behaviour today. Worth saying so in the PR description — this is a code/docs reconciliation, not a behaviour invention.

### CRITICAL — Step 3's `:699` repoint dereferences a process that was never invoked

The block gated at `:699` reads `BUILD_FULL_HYBRID_GTF.out.output` **directly**, twice:

- `main.nf:710` — `GFFREAD_FULL_HYBRID( BUILD_FULL_HYBRID_GTF.out.output, ch_fasta )`
- `main.nf:717` — `STAR_GENOMEGENERATE_FULL_HYBRID( ..., BUILD_FULL_HYBRID_GTF.out.output.map { ... } )`

Plan assumption 4 keeps `BUILD_FULL_HYBRID_GTF` gated on `extended_orf_active` (`:421`). So in a no-novel-source run the process is never invoked, and entering the `:699` block raises, at workflow-construction time:

```
Access to 'BUILD_FULL_HYBRID_GTF.out' is undefined since the process
'BUILD_FULL_HYBRID_GTF' has not been invoked before accessing the output attribute
```

This is reachable in the plan's headline configuration, and it is reachable **under `-profile test`**: `conf/test.config:28` sets `contrasts`, so `ch_contrasts_file` is non-null; `skip_plastid` defaults to false. Therefore:

- the plan's own new test (step 9) fails immediately;
- Validation 1 fails;
- Validation 4 (dotseq) fails — dotseq's ORF path lives inside this very block (`:788-799`);
- Validation 1's expectation "no `BUILD_FULL_HYBRID_GTF` tasks appear" cannot hold while `:710`/`:717` reference its output.

Step 4's instruction ("audit every use of `ch_full_hybrid_gtf`") **will not find these two lines** — they are not uses of `ch_full_hybrid_gtf`, they are direct process-output references. A grep for the channel name misses them. This is the single most likely way the implementation goes wrong.

Note the failure is loud, not silent, which limits the blast radius — but it means step 3 cannot land as three one-line substitutions.

Full audit of what the `:699` block needs (all four must be resolved, not two):

| Line | Reference | Needs |
|---|---|---|
| `:710` | `BUILD_FULL_HYBRID_GTF.out.output` | `[meta, gtf]`-shaped source |
| `:717` | `BUILD_FULL_HYBRID_GTF.out.output` | `[meta, gtf]`-shaped source |
| `:725` | `ch_full_hybrid_gtf` | value channel (per-sample STAR) |
| `:751` | `ch_full_hybrid_gtf` | GTF for Salmon tx2gene |

### CRITICAL — Step 5's dotseq fix is underspecified and, read literally, makes the error worse

`dotseqPrerequisitesError()` triggers on the **gate** at `utils_nfcore_riboseq_pipeline/main.nf:251`:

```groovy
if (!extended_orf_active || !any_caller_enabled || params.skip_plastid || !params.contrasts) {
```

The plan says only *"It currently lists 'a novel-transcript source' among dotseq's requirements … that clause must go"*. The clause is `:254` (`else if (!novel_source_configured) missing << ...`), which only builds the *message*. Deleting just that line leaves `:248`/`:251` intact, so with `--extended_orf_analysis true --skip_stringtie true` + callers + plastid + contrasts:

- `extended_orf_active` is still `false` → gate fires,
- every `if` inside is false → `missing == []`,
- the user gets `error("--translational_efficiency_method dotseq runs only at the ORF level and requires: . Pick anota2seq or deltate …")`.

A hard error with an empty requirement list is strictly worse than today's accurate-but-obsolete one. The fix must replace the predicate: drop `:247-248`, define `def orf_catalogue_active = params.extended_orf_analysis && any_caller_enabled`, gate on `if (!orf_catalogue_active || params.skip_plastid || !params.contrasts)`, and delete the `else if` at `:254` (turning `:253` back into a plain `if`). `novel_source_configured` then becomes unused in that function and must be removed too, or strict-syntax linting flags it.

### Important — the plastid warning's *guard* is the bug, not its wording

`utils_nfcore_riboseq_pipeline/main.nf:201`:

```groovy
if (params.extended_orf_analysis && novel_source_configured && params.skip_plastid) {
```

Plan step 7 asks whether the *wording* implies the catalogue is conditional on a novel source. It does not — `:202` already says "the ORF catalogue will still be built", which becomes *more* correct after the change. The real defect is the `novel_source_configured` conjunct: after the change, `--extended_orf_analysis true --skip_stringtie true --skip_plastid true` builds a catalogue with no ORF quantification and emits **no warning at all** — precisely the plan's own listed edge case, and precisely the silence that caused the original confusion. `docs/usage.md:664` promises that warning. Drop the conjunct.

### Important — the `ch_gtf` equivalence has one real exception: user-supplied `--canonical_gtf`

`--canonical_gtf` is a real user parameter (`nextflow.config:27`, `nextflow_schema.json:154`, threaded through `main.nf:94` → `prepare_genome/main.nf:131-136`). Two cases:

- **Derived** (`prepare_genome/main.nf:141-149`, AGAT longest-isoform + gffread): canonical transcript IDs are a subset of `ch_gtf`. Equivalence holds.
- **User-supplied** (`prepare_genome/main.nf:135`): an arbitrary GTF whose transcript IDs need not appear in `ch_gtf` at all.

Today's `full_hybrid` is immune to this because the awk appends hybrid rows whose `transcript_id` is absent from `ch_gtf`, and `ch_hybrid_gtf` ⊇ `ch_canonical_gtf` (`main.nf:311`) — so the catalogue GTF is **guaranteed** to contain every transcript ID any caller was given. A plain `ch_gtf` fallback loses that guarantee, and the consequence is silent, not loud. In non-extended mode Ribo-TISH and Ribotricer run against `ch_canonical_gtf` (`orf_caller_dispatch/main.nf:54-55, 135`), and the normaliser degrades on a missing ID rather than failing:

- `orfnormalise.py:557-563` (Ribo-TISH): *"ribotish predict reports one genomic span (GenomePos), not per-exon blocks, so the exon structure is recovered by intersecting that span with the transcript's exons"* — `if tx is None: blocks = [(start, end)]`. A multi-exon ORF collapses to one intron-spanning block.
- `orfnormalise.py:471-487` (RiboCode): falls back to `ORF_gstart`/`ORF_gstop`, else `continue` (row silently dropped).

An intron-spanning single-block BED12 then flows into `BEDTOOLS_GETFASTA -split -s -nameOnly` (`conf/modules.config:1092`) and `SEQKIT_TRANSLATE` → wrong nucleotide sequence, wrong peptide, wrong `aa_length`, wrong MMseqs clustering. So this is a genuine "silently produce wrong results" path. **The invariant "catalogue GTF ⊇ every enabled caller's annotation" is load-bearing and unstated in the plan.** Alternative A below satisfies it by construction; otherwise state it as an assumption and validate it.

### Important — no-RNA-seq runs silently produce no ORF-DTE

The RNA-seq-presence guard at `main.nf:316-322` fires only for `run_stringtie`. With no novel source, a Ribo-seq-only samplesheet plus `--contrasts` is perfectly legal. After step 3 the ORF-DTE block is constructed, `ch_rnaseq_reads` (`:720`) is empty, `FASTQ_ALIGN_STAR_FULL_HYBRID` and `QUANTIFY_HYBRID_RNA` yield nothing, `DTE_COUNTS_PREP` never runs — no output, no error. Behaviour item 5 ("ORF-level DTE runs when the catalogue exists, `--skip_plastid` false and `--contrasts` supplied") is therefore not accurate as stated. This partially pre-exists for `--novel_gtf` + Ribo-seq-only runs, but the change makes it far more reachable. Note that the Important-4 fix below turns it into a clear error: `dte_counts_prep.py:140-142` raises *"Secondary matrix has no sample columns left after dropping primary-role overlaps"*.

### Important — `orf_agreement_min_callers` is dead code; the plan's single-caller reasoning cites the wrong variable

`main.nf:463-465` (`orf_agreement_min_callers`), `:466` (`ch_enabled_orf_callers`) and `:467` (`ch_rank_aggregation_callers`) have **no consumers anywhere in the repo** — a full grep across `*.nf`/`*.config` returns only their definitions. The live consensus threshold is `params.orf_min_callers`, consumed by `CUSTOM_ORFMERGE` via `conf/modules.config:1081` and `:1104` (`--min-callers ${params.orf_min_callers} --min-samples ${params.orf_min_samples}`), default `1` (`nextflow.config:164-165`).

So the plan's Context bullet ("`orf_agreement_min_callers` at `:463` derives a strict majority … Both already handle a single caller") and its single-caller edge case describe arithmetic that never executes. The claim's *conclusion* is still fine — with `--orf_min_callers 1` the consensus view equals the full catalogue — but the reasoning should point at `params.orf_min_callers`. The real uncovered case: **`--orf_min_callers 2` with one enabled caller yields an empty consensus view, silently.** That combination becomes far more reachable now, and it is the default caller shape under `-profile test` (see below).

### Minor factual slips

- Self-Review: *"`extended_orf_active` used at seven sites in `main.nf` and again in two validators."* Seven occurrences in `main.nf` is right if the definition counts (`:362, 366, 421, 439, 477, 546, 699`). But `extended_orf_active` appears in only **one** validator (`:248`); it is `novel_source_configured` that appears in two (`:197`, `:247`). Immaterial to the design, but the second validator (`validateInputParameters`) needs a different edit from the first, so the distinction matters for step 5/6/7 sequencing.
- Open question about `#204` shifting `:699`: resolved — nothing moved.

---

## 2. Assumptions

| Assumption | Status |
|---|---|
| 1. Fallback is `ch_gtf`, not `ch_canonical_gtf` | **Confirmed** by `orf_caller_dispatch/main.nf:63-65` (same ternary already exists) and `assets/build_full_hybrid_gtf.awk`. One exception: user-supplied `--canonical_gtf` (Important above). |
| 2. Caller routing still needs a real novel source | **Confirmed** — `ch_hybrid_transcriptome_bam` is `channel.empty()` outside the second pass (`main.nf:364`), and RiboCode would get zero tasks. |
| 3. `--extended_orf_analysis` stays the catalogue switch | Holds, but see Alternative D on flag naming. |
| 4. `BUILD_FULL_HYBRID_GTF` stays gated; only the consumer gains a fallback | **This is the assumption that produces Critical 1.** It is not free: two consumers read the process output directly. Reconsider (Alternative A). |

Unstated assumptions the plan should surface:

1. **Catalogue GTF ⊇ every enabled caller's annotation.** Load-bearing (see Important above); currently guaranteed by the awk, not by `ch_gtf`.
2. **The ORF-DTE RNA denominator must be quantified against the catalogue GTF.** The comment at `main.nf:700-708` justifies the hybrid quant *solely* by novel genes lacking rows in the primary matrix. With no novel genes that justification is void, but the plan carries the machinery over unexamined.
3. **Channel type of `ch_gtf`.** `prepare_genome/main.nf:80` assigns a value channel (`Channel.value(file(gtf))`), but `:99` reassigns it from `CUSTOM_GTFFILTER.out.gtf.map{}` — a queue channel — whenever `filter_gtf` is true (`:93`, the default). `ch_full_hybrid_gtf` is explicitly `.first()`-ed (`main.nf:427`), so `ch_catalogue_gtf` would be a value channel on one branch and possibly a single-item queue on the other. This is **safe at `:487`** only because the catalogue subworkflow re-applies `.first()` (`orftable.../main.nf:44`). It is **not safe at `:725`** — `FASTQ_ALIGN_STAR_FULL_HYBRID` fans the GTF out per RNA-seq sample, and `orf_caller_dispatch/main.nf:37-38` documents exactly this contract: *"must be value channels … a single-item queue would serve only one sample."* If step 4's audit leads anyone to swap `:725` to a non-`.first()`-ed `ch_catalogue_gtf`, all but one sample silently vanish. Add `.first()` to the fallback branch.
4. **`params.extended_orf_analysis` is a real Boolean** — fine, nf-schema coerces, and `:362` already relies on it.

---

## 3. Efficiency analysis

**The plan's "no runtime cost … the fallback is a channel assignment" is wrong for any run with `--contrasts` and plastid enabled** — which is the default test profile and the common production shape.

After step 3, a no-novel-source run entering `:699` executes, against a GTF that is *identical to `ch_gtf`*:

| Process | Line | Cost | Already exists as |
|---|---|---|---|
| `GFFREAD_FULL_HYBRID` (`-w`) | `:709` | minutes | `ch_transcript_fasta` |
| `STAR_GENOMEGENERATE_FULL_HYBRID` | `:715` | **full genome index build** — tens of GB RAM, hours at human scale | `ch_star_index` |
| `FASTQ_ALIGN_STAR_FULL_HYBRID` | `:722` | second full STAR pass over every RNA-seq sample | `FASTQ_ALIGN_STAR` (`:242`) |
| `BAM_DEDUP_UMI_HYBRID_RNA` | `:734` | second dedup pass | `BAM_DEDUP_UMI` (`:268`) |
| `QUANTIFY_HYBRID_RNA` | `:746` | second Salmon quant + tximport | `QUANTIFY_STAR_SALMON` (`:566`) |

That dwarfs `MMSEQS_EASYCLUSTER`, which the Self-Review names as "the heaviest step". It also publishes `quantification/salmon_hybrid/` (`conf/modules.config:716, 726, 754, 762, 771`) with contents duplicating `quantification/salmon/` — user-visible confusion, and a much larger new snapshot.

**This is avoidable, and the codebase already anticipates it.** `modules/local/dte_counts_prep/templates/dte_counts_prep.py:17-21` and `:126-142` explicitly handle being handed *"a Salmon all-sample matrix carrying its Ribo-seq columns alongside the wanted RNA-seq ones"*: overlapping columns are dropped from the secondary matrix, leaving the RNA-seq columns. So when no novel source is configured, `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled` is a **drop-in** for `QUANTIFY_HYBRID_RNA.out.counts_gene_length_scaled` — same tximport shape, same annotation (`ch_gtf`, `main.nf:571`) — and the entire five-process chain can be skipped. It is also *strictly better* on the plan's own terms: the comment at `:704-708` wants canonical genes to "match the gene-level denominator", and reusing the primary matrix makes them match exactly rather than approximately.

Recommended: branch the denominator on `extended_orf_active` inside the `:699` block. That single branch fixes Critical 1 (the two direct `.out` references end up inside the novel-source arm) and removes the redundant compute in one move.

Everything else about the plan is cost-free as claimed: predicate evaluation and channel assignment happen at workflow construction.

---

## 4. Validation sufficiency

The six proposed validations are well-chosen but three cannot pass as written, and the highest-risk failure modes are uncovered.

**Cannot pass as written**

- **Validation 1** — crashes on Critical 1 under `-profile test` (contrasts set at `conf/test.config:28`). Its "no `BUILD_FULL_HYBRID_GTF` tasks appear" expectation is also design-dependent: under Alternative A, that process *should* appear (and publishes nothing without `--save_reference`, `conf/modules.config:718`, so it adds no snapshot noise).
- **Validation 3** — the Self-Review already concedes this is the weakest step. It is also unnecessary as prose: `tests/.nftignore_orf` does **not** ignore `orf_catalogue/**`, so `stable_path` already content-snapshots the catalogue. Assert catalogue *content* in the new test and diff the two `.snap` files — that is the concrete comparison the plan says it lacks, for free.
- **Validation 4** — crashes on Critical 1 for the same reason (dotseq's ORF path is inside the `:699` block at `:788-799`).

**Test-design corrections**

- `conf/test.config:33` sets **`skip_ribocode = true`**, so *every* pipeline test runs a **single caller (Ribo-TISH)**. This resolves the plan's internal inconsistency between step 9 ("at least one caller") and Validation 1 ("one caller") — under `-profile test` they are the same thing, and no caller-flag override is needed. It also means the single-caller edge case is already the only case any test has ever covered. Corollary: the comment in `tests/dotseq.nf.test:20-24` ("Default ribotish + ribocode is enough") is **stale** — ribocode is off in the test profile. Worth fixing while nearby.
- Step 9's assertion ("catalogue outputs exist") is too weak — an empty file passes. Assert non-zero size and let the snapshot lock content, mirroring `tests/stringtie_extended.nf.test:47-53`.
- No test registration is needed: `tests/tags.yml` maps `pipeline_riboseq` to `**.nf.test` and CI shards by discovery.

**Uncovered failure modes, in priority order**

1. **`BUILD_FULL_HYBRID_GTF.out` undefined** — Critical 1. Covered the moment the new test runs with contrasts, which it will by default. Do not weaken the new test by unsetting `contrasts` to dodge this.
2. **Catalogue GTF missing a caller's transcript IDs** (Important, user-supplied `--canonical_gtf`) — silent wrong peptides, no proposed check. Cheapest guard: assert in the new test that every `transcript_id` in `orf_catalogue/cohort.catalogue.tsv` resolves, or assert `orf_catalogue/normalised/*.bed12` block counts match the extended run for shared ORFs.
3. **Redundant hybrid quant chain** — Important 4. Assert `quantification/salmon_hybrid/` is **absent** in the new test. Cheap, and it locks the efficiency fix in place.
4. **Plastid warning suppression** — add a `--skip_plastid true` case asserting the warning fires (`workflow.stdout`/`workflow.trace` style assertion), since Validation 6 only checks the `:194` warning.
5. **`--orf_min_callers 2` with one caller** → empty consensus view, silent. Not worth a pipeline test; worth a `log.warn` in `validateInputParameters()`.
6. **Ribo-seq-only + contrasts** → no ORF-DTE, silent (Important 5). At minimum document the expected behaviour; better, let `dte_counts_prep.py:140-142` be the loud failure.

**Missing from Validation entirely:** `nf-core lint` (the schema text changes in Important 2 affect it) and `nextflow lint`/language-server strict-syntax checks (the unused `novel_source_configured` after step 5, and undeclared `ch_catalogue_gtf`).

---

## 5. Alternatives

### A. Run `BUILD_FULL_HYBRID_GTF` whenever the catalogue is active — no second GTF channel at all (recommended)

Change the guard at `main.nf:421` from `if (extended_orf_active)` to `if (extended_orf_active || orf_catalogue_active)` and delete plan step 2 and step 4 entirely.

With no novel source, `ch_hybrid_gtf == ch_canonical_gtf` (`main.nf:311`), so the awk emits `ch_gtf` plus any canonical row whose transcript is absent from `ch_gtf` — byte-identical to `ch_gtf` in the derived-canonical case, and *correct* in the user-supplied-canonical case.

Why this is better than the plan's design:

- **Fixes Critical 1 for free.** `BUILD_FULL_HYBRID_GTF.out.output` at `:710`/`:717` is always defined; the `:699` block needs no edit at all.
- **Satisfies the catalogue-GTF-superset invariant by construction**, closing the Important `--canonical_gtf` hole.
- **One code path.** No `ch_catalogue_gtf`, no channel-type asymmetry, no "which GTF did this run use?" ambiguity, no `.first()` question.
- **Smaller diff**: steps 2, 3 (the `:699` half) and 4 collapse into one guard change plus two `if`-condition swaps.

Costs: one extra awk task (a two-file streaming pass over a GTF — seconds), and a comment update at `orf_caller_dispatch/main.nf:43-44` ("empty and unread otherwise" becomes conditional). **No snapshot impact** — `conf/modules.config:718` publishes the output only under `--save_reference`. The plan's Validation 1 expectation about absent `BUILD_FULL_HYBRID_GTF` tasks needs rewording.

This does *not* address the redundant hybrid RNA quant; pair it with B.

### B. Branch the ORF-DTE RNA denominator on `extended_orf_active` (needed under either design)

Inside `:699`, when no novel source is configured, use `QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled` as `DTE_COUNTS_PREP`'s secondary matrix and skip `GFFREAD_FULL_HYBRID`, `STAR_GENOMEGENERATE_FULL_HYBRID`, `FASTQ_ALIGN_STAR_FULL_HYBRID`, `BAM_DEDUP_UMI_HYBRID_RNA` and `QUANTIFY_HYBRID_RNA`. Justified by `dte_counts_prep.py:17-21, 126-142`; removes hours of duplicated compute and a duplicated output directory; turns the empty-RNA case into a clear error. One caveat to state: the primary matrix comes from `ch_transcript_fasta` rather than a gffread extraction of the same GTF, so the two paths are equivalent at the annotation level but not necessarily bit-identical — which is fine, since the comparison baseline is "no output at all".

### C. The plan as written, plus the two Critical fixes

Viable. Keeps `ch_full_hybrid_gtf` strictly meaning "extended mode only", which has documentation value. But it needs a `ch_catalogue_gtf` **and** a resolution for `:710`/`:717` **and** an explicit statement of the superset invariant — strictly more moving parts than A for the same outcome.

### D. Worth raising, out of scope

The plan's step 1 says *"Name it for what it gates, not for the flag."* The same argument applies one level up: after this change `--extended_orf_analysis` means two unrelated things — "route novel transcripts into the callers" *and* "build the cross-caller catalogue". That is the same conflation being fixed internally, left in place at the user interface. A dedicated `--orf_catalogue` flag defaulting to `extended_orf_analysis`'s value would resolve the plan's first open question cleanly without changing default outputs. Not for this PR; worth capturing.

---

## 6. Action items

### Critical

1. **Resolve `BUILD_FULL_HYBRID_GTF.out.output` at `main.nf:710` and `:717` before repointing the `:699` gate.** Preferred: Alternative A (widen the `:421` guard). Otherwise the `:699` block must branch. Note explicitly in the plan that a grep for `ch_full_hybrid_gtf` does not find these two lines.
2. **Rewrite step 5 to change the *gate* at `utils_nfcore_riboseq_pipeline/main.nf:251`, not just the message at `:254`.** Deleting `:254` alone yields `error("… requires: .")`. Replace `:247-248` with `def orf_catalogue_active = params.extended_orf_analysis && any_caller_enabled`, gate on that, delete the `else if` at `:254`, and remove the now-unused `novel_source_configured`.

### Important

3. **Drop the `novel_source_configured` conjunct from the plastid warning guard at `utils_nfcore_riboseq_pipeline/main.nf:201`** — the wording is already correct; the guard silently suppresses the warning in the newly-enabled configuration that `docs/usage.md:664` promises it for. Restate step 7 accordingly.
4. **Add `nextflow_schema.json:806-807` to step 8.** Both the `description` ("Requires `--skip_stringtie false` or `--novel_gtf <path>`") and the `help_text` ("When true but no novel-transcript source is configured, the flag is a no-op") become wrong. Also update the in-code comment at `main.nf:870` ("empty unless extended-ORF + plastid both active"). Note that `docs/output.md:493, 504` and `docs/usage.md:431, 645, 654` already describe the post-change behaviour and need no edit — only `docs/usage.md:625` and the schema are actively wrong.
5. **Adopt Alternative B** (reuse the primary Salmon matrix as the ORF-DTE denominator when no novel source is configured) and correct the plan's Efficiency section, which currently claims "no runtime cost" for a change that can add a full STAR genome index build plus a second alignment and quantification pass.
6. **State the catalogue-GTF-superset invariant as an explicit assumption**, and either satisfy it by construction (Alternative A) or document that `--canonical_gtf` + `ch_gtf` fallback can silently collapse multi-exon ORFs to single blocks via `orfnormalise.py:557-563` / `:471-487`.
7. **Fix the test plan (step 9):** `conf/test.config:33` sets `skip_ribocode = true`, so `-profile test` is already single-caller — no caller override needed, and the step 9 / Validation 1 wording ("at least one" vs "one") should be unified. Assert non-zero catalogue sizes and let the content snapshot stand in for Validation 3. Add an assertion that `quantification/salmon_hybrid/` is **absent**. Fix the stale comment at `tests/dotseq.nf.test:20-24`.
8. **Correct the Context bullet about `orf_agreement_min_callers`.** `main.nf:463-467` is dead code with no consumers; the live threshold is `params.orf_min_callers` via `conf/modules.config:1081, 1104`. Add the uncovered case: `--orf_min_callers 2` with one enabled caller silently yields an empty consensus view — worth a `log.warn` in `validateInputParameters()`.

### Optional

9. Add `.first()` to the fallback branch so `ch_catalogue_gtf` is a value channel on both arms (if Alternative A is not taken). Safe at `:487` only because `orftable.../main.nf:44` re-applies it; unsafe anywhere per-sample, per `orf_caller_dispatch/main.nf:37-38`.
10. Declare `orf_catalogue_active` and `ch_catalogue_gtf` with `def` for strict-syntax cleanliness. Consider `as boolean` on `orf_catalogue_active` — it currently evaluates to a `List`, which is fine for `if` but not if it is ever passed as a subworkflow `val` (as `extended_orf_active` is at `:439`).
11. Reuse the existing idiom at `orf_caller_dispatch/main.nf:63-65` rather than introducing a second spelling of the same ternary.
12. Add a `CHANGELOG.md` entry — the repo convention is a `[#NNN]` bullet per behaviour change.
13. Note in "Remaining risks" that the predicate now has **two** definitions to keep in sync across three files (`main.nf`, and both validator functions), doubling the drift surface the plan already flags. Deferring the shared helper is reasonable; the cost has grown.
14. Add `nf-core lint` and `nextflow lint` to Validation — the schema edit and the unused-variable removal both touch what they check.
15. Correct the stale open question about `#204`: `HEAD = 6b17ac7` and every line number in the plan still resolves correctly.
