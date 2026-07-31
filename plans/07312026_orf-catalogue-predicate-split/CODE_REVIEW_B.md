# Code review B — `feat/orf-catalogue-without-novel-source` (commit `9f4a929`)

Reviewer B, fresh context. Repo `/Users/fkrueger/Github/riboseq`, diff `origin/dev..HEAD`
(single commit). Plan: `plans/07312026_orf-catalogue-predicate-split/PLAN.md`.

## Summary

The core change is correct and does what it claims. The predicate split is sound, the
widened `BUILD_FULL_HYBRID_GTF` guard genuinely is required rather than merely tidy, the
dead-code removal is verified safe, and the newly-reachable configuration constructs and
gates cleanly. `docs/output.md:493,504` already described the new gating, so this commit
closes a real doc/code divergence rather than inventing behaviour.

Two things need attention before merge:

1. **The new test ships without a snapshot.** CI runs `nf-test test --ci`
   (`.github/actions/nf-test/action.yml:58-68`), which fails rather than creates one.
   This is a guaranteed red build.
2. **The RNA-denominator substitution quietly drops the catalogue↔denominator gene-id
   namespace guarantee** that `docs/usage.md:654` asserts and that the commit message
   spends two paragraphs defending upstream. The plan anticipated the `--canonical_gtf`
   hazard for the *catalogue* (Validation 5) and fixed it there, but neither plan nor
   implementation followed it downstream into step 5. This is the shared blind spot.

Everything else is smaller: three duplicated encodings of the enabled-caller set, a few
plan items silently dropped, and some stale prose.

I applied three low-risk fixes (§Fixes applied) and left the rest as recommendations.

---

## Logic and scientific correctness

### 1. Superset invariant — holds in the default no-novel-source case, verified

`assets/build_full_hybrid_gtf.awk:13-26` prints file 1 (`ch_gtf`) verbatim including
comments, then appends only file-2 rows whose `transcript_id` (or, for transcript-less
rows, `gene_id`) is absent from file 1. With no novel source `ch_hybrid_gtf ==
ch_canonical_gtf` (`workflows/riboseq/main.nf:311`), so the output is `ch_gtf` plus
canonical-only transcripts.

Against the annotations each caller actually receives in that configuration:

| Caller | Annotation | Source |
| --- | --- | --- |
| Ribo-TISH | `ch_canonical_gtf` | `orf_caller_dispatch/main.nf:54-55,94-96` |
| Ribotricer | `ch_canonical_gtf` | `orf_caller_dispatch/main.nf:133-135` |
| Rp-Bp, PRICE | `ch_gtf` | `orf_caller_dispatch/main.nf:63-65` |
| RiboCode | `ch_gtf` | `orf_caller_dispatch/main.nf:221-223` |

`full_hybrid ⊇ ch_gtf ∪ canonical`, so the invariant holds. It is load-bearing: for
RiboCode's transcript-coordinate ORFs, a missing transcript falls through to
`orfnormalise.py:481-491`, which emits a **single genomic block** spanning
`ORF_gstart..ORF_gstop` — the silent intron-spanning degradation the commit message
describes. That claim checks out.

With `--canonical_gtf` unsupplied (the common case) canonical is AGAT+gffread-derived from
`ch_gtf` (`prepare_genome/main.nf:138-149`), so no rows are appended and `full_hybrid` is
byte-equivalent to `ch_gtf`. `BUILD_FULL_HYBRID_GTF` is then one gawk pass producing a
copy — accepted cost, correctly reasoned in the commit message.

### 2. Rp-Bp / PRICE routing is not disturbed — verified

`ch_full_hybrid_gtf` is now populated in configurations where it used to be
`channel.empty()`, but `orf_caller_dispatch/main.nf:63-65` selects
`ch_fasta_gtf_multi_isoform` on `extended_orf_active`, not on whether the GTF exists. So
Rp-Bp and PRICE still get `ch_gtf` without a novel source. No behaviour change.

### 3. Workflow construction is safe — verified

`BUILD_FULL_HYBRID_GTF.out.output` is read at `main.nf:692,699`, both inside
`if (extended_orf_active)` (`:690`), which implies the invocation at `:431-436`. Constructs
in every configuration.

`QUANTIFY_STAR_SALMON` is unconditional (`main.nf:551`), so the `.out` read in the else
branch (`:748`) always constructs. Its `counts_gene_length_scaled` now has three consumers
(`:573`, `:627`, `:748`); multiple consumers of a DSL2 process output are fine.

### 4. HIGH — the denominator substitution breaks the documented gene-id namespace guarantee

`docs/usage.md:654` (pre-change) asserted: *"The union is also the reference the ORF-level
DTE RNA denominator is quantified against, so catalogue gene ids share a namespace with
that matrix."* After this commit that is only true when a novel source is configured.

Without one:

- the catalogue resolves host genes against `full_hybrid` = `ch_gtf` ∪ canonical-only
  (`main.nf:472`);
- the denominator is the primary Salmon matrix, whose `tx2gene` derives from `ch_gtf`
  alone (`main.nf:551-565`, `:748`);
- `dte_counts_prep.py:150-156` keeps an ORF only if some candidate host gene is in
  `secondary.index`, and `:173-182` reports the drop count on stderr and exits 0.

So any gene id contributed by the canonical GTF but absent from `ch_gtf` puts its ORFs
outside the denominator's namespace and they are **silently dropped from ORF-level DTE**.
The commit message cites exactly this configuration (user-supplied `--canonical_gtf`,
"MANE Select, whose IDs are absent from an Ensembl GTF") as the reason to build the full
hybrid GTF unconditionally — and then reintroduces the mismatch one stage downstream.

Blast radius, honestly bounded:

- **Total** namespace divergence is already broken elsewhere. `FILTER_COUNTS_CANONICAL`
  keys the gene-level matrix on canonical `gene_id`s
  (`assets/filter_counts_canonical.awk:13-22`), so a fully foreign canonical GTF already
  yields a header-only gene-level TE matrix. A user in that state notices.
- **Partial** divergence is the live case and is not obviously broken: a `--canonical_gtf`
  from a different Ensembl release than `--gtf`, or a MANE GTF with version-suffixed ids
  (`ENSG…​.5` vs Ensembl GTF's unsuffixed `gene_id` + separate `gene_version`). Gene-level
  TE still mostly works, `full_hybrid` gains the extra genes, and the ORFs Ribo-TISH and
  Ribotricer call on them vanish from ORF-level DTE with only a stderr line. Wrong-but-
  plausible: an ORF-level DTE table silently missing one subset of callers' ORFs.

`tests/stringtie_extended.nf.test:34-43,51-52` already asserts precisely this invariant
for the novel-source path (`catalogue_genes.every { rna_genes.contains(it) }`). The new
test has no equivalent. Recommendation and exact code in §R1.

### 5. MEDIUM — TI-seq columns, and an asymmetry between the two denominators

`tiseq` is a supported sample type (`assets/schema_input.json:38-39`). The two denominators
carry different sample sets:

- novel-source path: `QUANTIFY_HYBRID_RNA` is fed RNA-seq reads only (`main.nf:702-713`),
  so its matrix has **rnaseq columns only**;
- no-novel path: the primary Salmon matrix covers all sample types; `dte_counts_prep.py:129-137`
  drops only the columns overlapping primary (the Ribo-seq ones), leaving **rnaseq + tiseq**.

Consequences, both directions:

- The no-novel path is the more robust one. `deseq2_deltate.R:379-382` hard-stops when any
  samplesheet row is missing from the count table, so a samplesheet containing `tiseq`
  samples makes `DESEQ2_DELTATE_ORF` **fail** on the novel-source path today. That is a
  pre-existing bug in #204, not introduced here, and no test covers it (the chr20 test data
  has no tiseq samples).
- Conversely, ORF-level DTE results are not comparable between the two configurations
  beyond the gene set, because the extra columns enter the deltaTE interaction model as a
  third `type` factor level (`deseq2_deltate.R:396-405,414-428`).

Not a blocker for this commit; worth a line in the docs and an issue for the novel path.

### 6. Ribo-seq-only + `--contrasts` turns a silent skip into a hard failure

Previously the novel path produced an empty `QUANTIFY_HYBRID_RNA` and `DTE_COUNTS_PREP`
never ran. Now `dte_counts_prep.py:138-143` raises `SystemExit` with a clear message. This
is an improvement (and the gene-level path already fails at `deseq2_deltate.R:400-402` in
that configuration), so no action — noting it because it is a behaviour change the commit
message does not mention.

### 7. tximport length scaling — no material issue

`counts_gene_length_scaled` is `lengthScaledTPM`, whose per-gene mean effective length is
averaged over every sample in the tximport call. The all-sample primary matrix therefore
does not produce bit-identical RNA-seq numbers to an RNA-only tximport, even for identical
genes. But the factor is per-gene and constant across samples, so it is absorbed by the
gene-wise intercept and does not bias contrasts. The commit's stated goal — the ORF-level
denominator becomes *exactly* the gene-level one — is achieved. No action.

---

## Errors and validation

### 8. CRITICAL — `tests/orf_catalogue_no_novel_source.nf.test` has no snapshot

The test calls `snapshot(...).match()` (`:50-57`) but no
`tests/orf_catalogue_no_novel_source.nf.test.snap` exists. CI runs:

```
nf-test test --profile=+docker --ci --changed-since HEAD^ ...
```

(`.github/actions/nf-test/action.yml:58-68`). `--ci` fails the run when a snapshot has to
be created. Every other enabled test in `tests/` ships its `.snap`. Generate it with a
local `-profile test,docker` run and commit it.

### 9. MEDIUM — the `!salmon_hybrid.exists()` assertion passed for the wrong reason

`BUILD_FULL_HYBRID_GTF` publishes into `${outdir}/quantification/salmon_hybrid`
(`conf/modules.config:712-720`) and is now *expected* to run in this configuration. The
original assertion survived only because `saveAs` returns `null` unless
`params.save_reference` (`:718`, default `false` at `nextflow.config:26`), so the directory
is never created. It would have failed under `--save_reference true` and could not
distinguish "hybrid quant skipped" from "reference not saved".

A concurrent reviewer in this session has already tightened this in the working tree, to
assert on `quantification/salmon_hybrid/salmon.merged.hybrid.gene_counts_length_scaled.tsv`
and `.../star_index` instead. I did not duplicate the change; the file names match
`tests/stringtie_extended.nf.test.snap`. Independently reached, so treat as confirmed.

Separately: `BUILD_FULL_HYBRID_GTF`'s publish path is now misleading — the GTF is a
catalogue input, and in this configuration nothing Salmon-hybrid happens. Consider moving
it to `orf_catalogue/` or `reference/`. Renaming the *process* would move every pipeline
snapshot (`nf_core_riboseq_software_mqc_versions.yml`), so change only the path.

### 10. MEDIUM — `--extended_orf_analysis true` with zero callers is still silent

`orf_catalogue_active` is false when no caller is enabled, so no catalogue, no
`orf_quantification/` and no ORF-level DTE — with no warning. The plan lists this edge case
("**No callers enabled** — no catalogue. Unchanged.") but "unchanged" is exactly the silent
no-op this commit exists to eliminate, and the sibling `--orf_min_callers` warning covers
the adjacent case. **Fixed** (§Fixes applied).

---

## Structure

### 11. MEDIUM — the enabled-caller set now has three encodings in two files

| Location | Shape |
| --- | --- |
| `workflows/riboseq/main.nf:361-366` | `List<String>` |
| `utils_nfcore_riboseq_pipeline/main.nf:199-205` | count of booleans |
| `utils_nfcore_riboseq_pipeline/main.nf:254` | boolean OR chain |

The commit added the second. Adding a sixth caller now needs three edits in three shapes,
and nothing fails if one is missed. The plan's Self-Review accepts the duplication because
validators run before the workflow and cannot import a workflow-scope `def` — true, but a
plain `def enabledOrfCallers()` helper in `utils_nfcore_riboseq_pipeline/main.nf`
(file-scope, alongside `validateInputParameters()`) is importable from
`workflows/riboseq/main.nf` and would collapse all three to one list plus `.size()`. Worth
doing; not blocking.

### 12. LOW — misc

- `main.nf:372` — `as boolean` is a no-op: Groovy's `&&` already yields a `boolean`, so
  `params.extended_orf_analysis && enabled_orf_callers` is never a `List`. The plan's
  justification for the cast (step 1: *"evaluates to a `List`"*) is incorrect, and the
  inconsistency with `:371`, which needs no cast for the same reason, invites the question.
  Harmless; drop it or cast both.
- `main.nf:344-366` — the caller enumeration is wedged between the "Extended ORF discovery:
  second STAR pass" comment block and the `if (extended_orf_active)` block at `:376` that
  the comment describes. Move it above the comment.
- `orf_caller_dispatch/main.nf:43` — the new comment drops the useful half of the old one.
  `ch_full_hybrid_gtf` is now populated whenever the GTF is built but still *read* only
  when `extended_orf_active` (`:63-65`). Suggest: `// populated whenever the hybrid GTF is
  built; read only in extended mode.`
- `main.nf:859` — `orf_count_matrix … empty unless extended-ORF + plastid both active` is
  the stale emit comment the plan flagged (as `:870`) and it was not updated. The gate is
  now `--extended_orf_analysis` + a caller + plastid.
- `CHANGELOG.md:46` sits between `#204` and `#210`, out of numeric order. Cosmetic.

### 13. LOW — commit-message claim that does not hold

The commit message says `dotseqPrerequisitesError` "would have errored with an empty
requirements list". It could not: the old gate fired only if one of four conditions held,
and each of the four appended an item — including the `if/else if` pair at the old `:253-254`,
which covered both ways `extended_orf_active` could be false. The real fix is the
behavioural one (dotseq no longer demands a novel source), which is correct and worth
stating on its own. The plan's step 6 carries the same wrong premise; the resulting code is
right either way, and `missing` still cannot be empty (`main.nf:257-263`).

---

## Plan coverage

Steps 1-9 and 11 landed. Verified dropped:

- **Step 10, `docs/usage.md:439`** — listed explicitly, conditional on step 5 landing, which
  it did. Not updated: the page still stated unconditionally that the denominator is the
  novel-augmented transcriptome published under `quantification/salmon_hybrid/`. **Fixed.**
- **Step 10, cross-caller asymmetry note** — the plan asked for a doc note that Ribo-TISH and
  Ribotricer see the one-transcript-per-gene canonical GTF while Rp-Bp/PRICE/RiboCode see
  multi-isoform `ch_gtf`, so an isoform-specific ORF can never reach 2-of-2 agreement. The
  paragraph added at `docs/usage.md:627` covers `orf_id` instability and the single-caller
  consensus but not this. `:654` describes the asymmetry without its `--orf_min_callers`
  consequence. Not fixed — see §R2.
- **Step 10, `main.nf:870`** stale emit comment — not updated (§12).
- **Step 11, `tests/dotseq.nf.test:18-24`** stale comment — not updated; it still says
  "Default ribotish + ribocode is enough", but `conf/test.config:33` sets
  `skip_ribocode = true`.
- **Validation 3 and 5** — neither became a test. V5 (`--canonical_gtf` with ids absent from
  `--gtf`) is the one that matters: the plan flagged it as fiddly and warned that skipping it
  leaves the superset invariant resting on reading the awk. It was skipped, and §4 is the
  cost.

Nothing was implemented that the plan did not ask for.

---

## Fixes applied

Uncommitted, in the working tree. Not committed or pushed.

1. **`docs/usage.md:439`** — split into three paragraphs: which matrix supplies the
   denominator now depends on the novel source. Documents that without one the primary
   Salmon matrix (`quantification/salmon/`) is reused, that the ORF-level denominator is
   then exactly the gene-level one, and that the join drops the Ribo-seq columns it carries.
   (Missed plan item.)
2. **`docs/usage.md:654`** — scoped the namespace claim to the novel-source path and added
   the `--canonical_gtf` caveat from §4, including that the drop count surfaces on
   `DTE_COUNTS_PREP`'s stderr. Documentation only; §R1 is the real remedy.
3. **`utils_nfcore_riboseq_pipeline/main.nf:199-209`** — hoisted `enabled_caller_count`
   above both warnings and added one for `--extended_orf_analysis true` with every caller
   disabled (§10). Six added lines, matching the sibling warnings' style.

**Working-tree note.** A local `PostToolUse` formatter in my session reformatted all 175
lines of `utils_nfcore_riboseq_pipeline/main.nf` on my first edit (`log.warn "…"` →
`log.warn(…)`, `return` → `return null`, collapsed multi-line string literals). I reverted
it and re-applied only my six lines via a non-hooked path; `git diff` on that file is now
exactly the intended hunk. This is not a repo problem — `.pre-commit-config.yaml` runs
`nextflow-lint` in check mode only (`-output json`, no `-format`) — but be aware that a
`nextflow lint -format` run over this file will produce a ~170-line unrelated diff.

Also present in the working tree: another reviewer's tightening of the test assertions
described in §9.

---

## Recommendations

### R1 — Critical / High

1. **Commit a snapshot for the new test** (§8). Blocking.
2. **Add the namespace assertion to `tests/orf_catalogue_no_novel_source.nf.test`**, mirroring
   `tests/stringtie_extended.nf.test:34-43,51-52`. This is the guard for §4 and it is cheap,
   since the run already produces both files:

   ```groovy
   // Every catalogue host gene must have a row in the reused primary Salmon matrix,
   // or its ORFs are silently dropped from ORF-level DTE.
   def orf_to_gene = file("${params.outdir}/orf_catalogue/cohort.catalogue.orf_to_gene.tsv").readLines()
   def catalogue_genes = orf_to_gene.drop(1)
       .collect { it.split('\t', -1)[1] }
       .findAll { it }
       .unique()
   def rna_genes = file("${params.outdir}/quantification/salmon/salmon.merged.gene_counts_length_scaled.tsv")
       .readLines()
       .drop(1)
       .collect { it.split('\t')[0] }
       .toSet()
   ...
   { assert catalogue_genes.size() > 0 },
   { assert catalogue_genes.every { rna_genes.contains(it) } },
   ```

   I did not apply this: it cannot be verified without running the pipeline, and if it fails
   it has found a real gap rather than a test bug. Run it before deciding.
3. **Decide how to handle a divergent `--canonical_gtf`** (§4). Options, in ascending cost:
   (a) documentation only — applied; (b) `log.warn` when `--canonical_gtf` is user-supplied
   and no novel source is configured, since that is the only way the namespaces can diverge;
   (c) fall back to the `QUANTIFY_HYBRID_RNA` path when `params.canonical_gtf` is set, which
   restores the guarantee by construction at the cost of the STAR index build the commit set
   out to avoid — i.e. gate step 5's reuse on `!params.canonical_gtf` rather than on
   `extended_orf_active`. (b) is the proportionate choice.

### R2 — Medium

4. Extract a single `enabledOrfCallers()` helper and collapse the three encodings (§11).
5. Add the cross-caller asymmetry doc note the plan asked for, tied to `--orf_min_callers`
   (§Plan coverage).
6. Note the TI-seq asymmetry in the docs and open an issue for the novel-source path's
   `DESEQ2_DELTATE_ORF` failure on tiseq samplesheets (§5).
7. Move `BUILD_FULL_HYBRID_GTF`'s publish path out of `quantification/salmon_hybrid/`
   (path only, not the process name) (§9).
8. Say in the CHANGELOG that existing `--extended_orf_analysis true` users get new outputs
   *and* new compute on upgrade — the full catalogue subworkflow including
   `MMSEQS_EASYCLUSTER`, `QUANTIFY_ORF_PSITE` per Ribo-seq sample, and ORF-level DTE. It is
   the intended fix and gated behind a non-default flag, but it is a surprise for anyone
   already passing the flag.

### R3 — Low

9. `main.nf:859` stale emit comment; `tests/dotseq.nf.test:18-24` stale comment;
   `orf_caller_dispatch/main.nf:43` comment; `as boolean` at `main.nf:372`; move the caller
   enumeration above the second-STAR-pass comment; `CHANGELOG.md:46` ordering (§12).
10. Drop the "empty requirements list" claim from the commit message (§13).
