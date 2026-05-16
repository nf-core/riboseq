# PRICE (#170) implementation notes

Drafted 2026-05-16 in worktree `riboseq-mod-170-price` while v13 Platform iteration was running. Branch: `feat/170-price` off `feat/modernisation-aggregation@25597f2`.

## Status

- Bioconda recipe `gedi 1.0.6a` (noarch, ~22 MB) merged 2026-05-16 (bioconda-recipes #65445)
- Wave container verified: `wave.seqera.io/wt/c85a0e0481b1/wave/build:gedi-price--ef135d2a52782eb4`. **MUST pin `openjdk=17`** in the conda env — Java 25 (current conda-forge HEAD) SIGSEGVs at class discovery
- Biocontainer auto-publish still propagating to quay.io / Galaxy depot at time of writing. Use Wave on Seqera for now.
- environment.yml minimum:
  ```yaml
  channels: [bioconda, conda-forge]
  dependencies:
    - bioconda::gedi=1.0.6a
    - conda-forge::openjdk=17
  ```

## CLI surface

PRICE is a multi-step pipeline (codon stats → model estimation → ORF inference → FDR correction → reassignment) all driven by one command.

```
gedi -e Price <options>
```

**Required**

- `-reads <reads>` — input ribo-seq reads (BAM or CIT)
- `-prefix <prefix>` — output prefix (most artefacts named `${prefix}.<suffix>`)
- `-genomic <genomic>` — GEDI genome index (`.oml` file, built by `IndexGenome`)

**Useful**

- `-nthreads <N>` — default 10
- `-fdr <f>` — FDR cutoff (default 0.1)
- `-filter <a:b>` — read-length filter (e.g. `28:30`)
- `-novelTranscripts` — infer novel transcripts from reads (extended-mode flag!)
- `-percond` — estimate models per condition
- `-plot` — R plots
- `-seed <int>` — EM seed (default 1337)
- `-skipmt`, `-keepAnno`, `-inferDelta`, `-opt` — various toggles
- `-progress` — progress output

## Primary output

`${prefix}.orfs.tsv` — table containing all ORFs (this is what `ORF_NORMALISE_PRICE` will consume).

Other artefacts (mostly intermediate but worth publishing for inspection): `.codons.cit`, `.codons.bin`, `.estimateData`, `.model`, `.clusters.cit`, `.start.model`, `.noise.model`, `.orfs.pvals`, `.orfs.bin`, `.orfs.cit`, `.signal.tsv`.

## Genome index step

```
gedi -e IndexGenome -s reference.fa -a reference.gtf -nomapping -n hybrid_reference
```

- `-nomapping` skips kallisto/bowtie/STAR indices (we have STAR already)
- `-n <name>` sets the genome name; the index goes to `${PWD}/${name}.oml` by default (or `~/.gedi/genomic/${name}.oml`) plus a sibling directory of binary files

Verify directory layout against an actual run before writing the module. Outputs land in/around the OML directory.

## Suggested workflow structure (mirrors Rp-Bp pattern)

```
PRICE_INDEXGENOME (once per pipeline)
  ↓ oml + indices
PRICE_PRICE (one call, all sample BAMs together — see "samples handling" below)
  ↓ ${prefix}.orfs.tsv
ORF_NORMALISE_PRICE (per-sample? or once?)
  → joins build_orf_catalogue alongside ORF_NORMALISE_{RIBOTISH,RIBOCODE,RIBOTRICER,RPBP}
```

## Samples handling

PRICE expects all samples passed at once (multi-condition design like `gedi -e Price -reads samples.bamlist`). Unlike ribotish_predict which can be per-sample, PRICE estimates a shared codon-position model across samples. Wiring should be:

- gather all riboseq BAMs (the ones already going to `RIBOTISH_PREDICT_ALL`)
- emit a `samples.bamlist` listing them
- call PRICE_PRICE once

For `ORF_NORMALISE_PRICE` we may want to fan out the resulting `.orfs.tsv` per sample using PRICE's per-sample columns (if it has them) or just have one normalisation run that emits BED12 for the merged ORF set.

## Open questions for implementation

1. **Output column layout of `${prefix}.orfs.tsv`** — need to actually run PRICE on a small test BAM to see. Wiki mentions the columns but not in machine-readable form. Plan: run PRICE on chr20 test data in the next session and inspect the TSV header + first few rows.
2. **Multi-sample columns** — does `.orfs.tsv` carry per-sample counts/scores, or just merged? Affects whether `ORF_NORMALISE_PRICE` fans out or not.
3. **Memory profile** — PRICE is Java EM heavy. The wiki suggests cluster usage. Set a `label 'process_high_memory'` or similar, and probably `task.memory` >= 32 GB. Confirm with the small test run.
4. **`-novelTranscripts`** — wire under `extended_orf_analysis` so novel ORF discovery kicks in when StringTie is feeding hybrid GTF. Match the pattern used for ribotish/ribotricer.
5. **GEDI index reusability** — the `.oml` index could conceivably be cached and reused across runs of the same reference. For now, treat it as a per-pipeline output (similar to `STAR_GENOMEGENERATE`).

## Java version pin gotcha

Bioconda recipe declares `openjdk >=11` which lets conda solve to the newest, currently Java 25. Java 25 SIGSEGVs gedi at class-discovery time. Java 17 (LTS) works. Pin in environment.yml. **Consider upstreaming a `<22` upper bound to the Bioconda recipe** once we confirm which JDK versions are actually compatible.

## Version capture

`gedi -e Price -h` prints `Price version 1.0.4` — that's PRICE's internal version, not the package. Use `gedi -e Version` (returns `Gedi 1.0.6a`) for `versions.yml`. Don't grep the Price banner.

## Stub vs real module strategy

For module-level nf-test, the stub block should `touch` all expected outputs. For the full test, use a chr20 sliced BAM (from nf-core/test-datasets/riboseq) and confirm `.orfs.tsv` is non-empty.

## Next-session actionable list

1. Build a minimal local test: chr20 FASTA + chr20 GTF + one chr20 BAM. Run IndexGenome then Price end-to-end in Docker. Capture `${prefix}.orfs.tsv` header + first 10 rows.
2. Write `modules/local/price/indexgenome/main.nf`.
3. Write `modules/local/price/price/main.nf`.
4. Write `bin/orf_normalise_price.py` based on the captured columns.
5. Wire into `workflows/riboseq/main.nf` under `if (params.run_price)` (default false) and the `BUILD_ORF_CATALOGUE` subworkflow.
6. Add module + pipeline-level tests.
7. Update docs (`usage.md`, `output.md`, `README.md`, `CHANGELOG.md`).

## Working files captured here

- `.wave-build/gedi.yml` — environment.yml for the Wave container
- This notes file

No source code committed yet — only this planning doc. The next session can pick up from here cold.
