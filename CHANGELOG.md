# nf-core/riboseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev

### `Added`

- [#170](https://github.com/nf-core/riboseq/issues/170) - Add PRICE (Erhard et al. 2018) as an opt-in ORF caller via `--run_price` (default `false`). Invoked one-shot across the riboseq cohort (`gedi -e IndexGenome` + `gedi -e Price`); calls flow into the cross-caller ORF catalogue. Container via `bioconda::gedi=1.0.6a` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#125](https://github.com/nf-core/riboseq/pull/125) - Add rRNA removal tool selection with support for SortMeRNA (default), Bowtie2, and RiboDetector ([@pinin4fjords](https://github.com/pinin4fjords))
- [#128](https://github.com/nf-core/riboseq/pull/128) - Add DESeq2-based deltaTE analysis as an alternative to anota2seq for translational efficiency analysis ([@pinin4fjords](https://github.com/pinin4fjords))
- [#131](https://github.com/nf-core/riboseq/pull/131) - Add ribotish quality output routing to MultiQC ([@pinin4fjords](https://github.com/pinin4fjords))
- [#135](https://github.com/nf-core/riboseq/pull/135) - Add optional read length equalisation to trim RNA-seq reads to match Ribo-seq lengths before quantification ([@pinin4fjords](https://github.com/pinin4fjords))
- [#136](https://github.com/nf-core/riboseq/pull/136) - Add RiboCode ORF detection with P-site analysis and metaplots ([@JackCurragh](https://github.com/JackCurragh))
- [#137](https://github.com/nf-core/riboseq/pull/137) - Add Salmon pseudo-alignment pathway for translational efficiency analysis ([@pinin4fjords](https://github.com/pinin4fjords))
- [#138](https://github.com/nf-core/riboseq/pull/138) - Add MultiQC configuration for BBSplit, Bowtie2 rRNA removal, UMItools, and UMIcollapse ([@pinin4fjords](https://github.com/pinin4fjords))
- [#140](https://github.com/nf-core/riboseq/pull/140) - Add P-site identification using plastid ([@suhrig](https://github.com/suhrig))
- [#141](https://github.com/nf-core/riboseq/pull/141) - Add in-frame P-site quantification method ([@suhrig](https://github.com/suhrig))
- [#142](https://github.com/nf-core/riboseq/pull/142) - Add `ribodetector_chunk_size` parameter to control RiboDetector memory usage ([@pinin4fjords](https://github.com/pinin4fjords))
- [#144](https://github.com/nf-core/riboseq/pull/144) - Add bigWig coverage tracks ([@suhrig](https://github.com/suhrig))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Add riboWaltz QC plots to the MultiQC report (P-site region distribution, reading-frame distribution, start/stop-codon metaprofiles) via the new MultiQC riboWaltz module ([MultiQC#3465](https://github.com/MultiQC/MultiQC/pull/3465)) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#157](https://github.com/nf-core/riboseq/issues/157) - Add optional StringTie reference-guided novel transcript discovery; merged GTF published as a side product (downstream wiring deferred to a follow-on, see PR #158) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#161](https://github.com/nf-core/riboseq/issues/161) - Add a one-transcript-per-gene canonical annotation backbone via `--canonical_gtf`, used by the genome-coordinate ORF callers, plastid P-site quantification and translational-efficiency analysis; falls back to AGAT longest-CDS isoform when not supplied ([@pinin4fjords](https://github.com/pinin4fjords))
- [#162](https://github.com/nf-core/riboseq/issues/162) - Document that Ribo-TISH `quality` continues to consume the canonical annotation (P-site offset calibration uses CDS-bearing canonical transcripts; mixing CDS-absent novel transcripts would degrade calibration without adding diagnostic signal) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Prefer RNA-seq BAMs for StringTie assembly when available; fall back to Ribo-seq BAMs with tightened parameters via the new `--stringtie_ribo_fallback_args` parameter and emit a runtime warning ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Add `--novel_gtf` to bypass StringTie and feed a user-supplied novel-transcript GTF directly into the filtering chain (gffcompare class-code filter + optional blacklist + hybrid concatenation) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Classify novel transcripts with gffcompare against the full reference and filter to user-configurable class codes via the new `--stringtie_class_codes` parameter (default `u`, intergenic only; stranded users may add `x` for antisense) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Add optional strand-aware rRNA/repeat blacklist intersect via the new `--rrna_blacklist` parameter (`bedtools intersect -v -s`) to remove novel-transcript assembly artefacts ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Build a hybrid GTF (`<outdir>/stringtie/hybrid_reference.gtf`) by concatenating the canonical backbone with the filtered novel transcripts; expose it as the `hybrid_gtf` workflow emit channel for downstream ORF-caller wiring ([@pinin4fjords](https://github.com/pinin4fjords))
- [#165](https://github.com/nf-core/riboseq/issues/165) - Add `--extended_orf_analysis` (default `false`) to wire the hybrid GTF into the genome-BAM ORF callers: Ribo-TISH `predict` receives the hybrid GTF on `-g` and the canonical backbone on `-a` for background and classification; Ribotricer `prepare-orfs` receives the hybrid GTF directly. RiboCode, riboWaltz and Salmon continue on the reference transcriptome due to the transcriptome-BAM architectural constraint, and plastid P-site quantification on the canonical backbone (follow-on second-STAR-pass for RiboCode tracked in #171). When the flag is enabled without a novel-transcript source (`--skip_stringtie false` or `--novel_gtf`), the pipeline warns and falls back to canonical ([@pinin4fjords](https://github.com/pinin4fjords))
- [#171](https://github.com/nf-core/riboseq/issues/171) - Under `--extended_orf_analysis true`, run a second STAR pass for Ribo-seq samples against a hybrid transcriptome (canonical + filtered novel intergenic) so RiboCode can call ORFs on novel transcripts. The new local `EXTENDED_ORF_SECOND_PASS_ALIGN` subworkflow extracts spliced transcript sequences from the hybrid GTF with `gffread -w` and rebuilds the STAR index with the hybrid GTF as `--sjdbGTFfile`; the hybrid transcriptome FASTA and index are built once per run and the second STAR pass runs only on Ribo-seq samples. RiboCode then consumes the hybrid transcriptome BAM and the hybrid GTF as its annotation source. riboWaltz and Salmon stay on the primary reference transcriptome and plastid P-site quantification on the canonical backbone by design. Hybrid alignment outputs land under `<outdir>/hybrid_star/`. The default-off path is unchanged ([@pinin4fjords](https://github.com/pinin4fjords))
- [#169](https://github.com/nf-core/riboseq/issues/169) - Add Rp-Bp as an opt-in Tier-1 ORF caller via `--run_rpbp` (default `false`). Rp-Bp is a Bayesian-strict caller (Malone et al. 2017) that complements RiboCode's permissive canonical-CDS calls; benchmark data (FK/NGB, May 2026) places it at Tier-1 for both rank concordance (mean Spearman 0.893) and set overlap (mean Jaccard 0.673). The pipeline pulls in the upstream `nf-core/rpbp/*` modules (`preparegenome`, `extractmetageneprofiles`, `estimatemetagenebayesfactors`, `selectperiodicoffsets`, `getperiodiclengthsoffsets`, `extractorfprofiles`, `estimateorfbayesfactors`, `selectfinalpredictionset`) wired into the `FASTA_GTF_BAM_RPBP` nf-core subworkflow. CLI overrides are exposed via `--extra_rpbp_preparegenome_args` and `--extra_rpbp_predictorfs_args`. Per-sample final-prediction BED + nucleotide/protein FASTAs land under `<outdir>/orf_predictions/rpbp/` and feed the cross-caller rank-aggregation set via the ORF catalogue's `CUSTOM_ORFNORMALISE`. The same conditional canonical-vs-hybrid GTF selection as Ribo-TISH/Ribotricer is applied. Expect ~20-24h per replicate at genome-wide scale; the pipeline emits a runtime warning when enabled ([@pinin4fjords](https://github.com/pinin4fjords))
- [#167](https://github.com/nf-core/riboseq/issues/167) - Add a cohort-level ORF catalogue under `--extended_orf_analysis true`. The `CUSTOM_ORFNORMALISE` upstream module converts each per-sample, per-caller output (Ribo-TISH, RiboCode, Ribotricer, Rp-Bp, PRICE) to a unified BED12 + sidecar TSV; the `ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE` nf-core subworkflow then merges them with a class-aware strategy (transcript-ID grouping for annotated multi-exon CDS, 80% reciprocal overlap for single-exon novel intergenic and smORFs <= 100 aa) and emits `*.catalogue.bed12`, `*.catalogue.tsv` (per-ORF table with `called_by_<caller>` / `score_<caller>` columns), `*.orf_to_gene.tsv`, `*.catalogue.aa.fasta` (AA sequences, via `bedtools getfasta` + `seqkit translate`) and a MultiQC custom-content per-class count table under `<outdir>/orf_catalogue/`. The catalogue runs once per pipeline invocation and gates on `--extended_orf_analysis true` plus a non-empty enabled-caller set; the default-off path is unchanged ([@pinin4fjords](https://github.com/pinin4fjords))
- [#166](https://github.com/nf-core/riboseq/issues/166) - Add per-ORF in-frame P-site quantification additive to the existing gene-level path. A new `QUANTIFY_ORF_PSITE` subworkflow expands the cohort BED12 catalogue into codon-start positions (frame defined by each ORF's own ATG, not GTF `phase`), runs per-sample `bedtools intersect` against the plastid wiggle tracks and assembles an ORF x sample count matrix (`orf_psite_counts.tsv`) zero-filled for ORFs absent from a sample. Runs whenever `--extended_orf_analysis true`, at least one ORF caller is enabled, and plastid is not skipped; warns and skips when `--skip_plastid true`. The matrix is published under `<outdir>/orf_quantification/` and emitted on the `orf_count_matrix` workflow channel for the future DTE step (#168). The gene-level `QUANTIFY_INFRAME_PSITE_PLASTID` path is unchanged ([@pinin4fjords](https://github.com/pinin4fjords))
- [#168](https://github.com/nf-core/riboseq/issues/168) - Add ORF-level differential translation analysis on top of the existing gene-level anota2seq / deltaTE DTE. Under `--extended_orf_analysis true` with `--te_quantification_method plastid_psite`, the gene-level TE Ribo-seq numerator is re-aggregated from per-ORF P-site counts summing ONLY `canonical_cds` ORFs (Tier 1), keeping uORF / dORF / novel_u / smORF dynamics out of the gene-level sum. A method-agnostic prep step expands the gene-level RNA-seq matrix onto ORF rows via `orf_to_gene.tsv`; the user's `--translational_efficiency_method` (anota2seq, deltate or dotseq) is then fit at ORF resolution (Tier 2), with the Ribo-seq numerator from `orf_psite_counts.tsv` and the RNA-seq denominator joined from gene-level Salmon counts (novel intergenic ORFs with no host gene are dropped, low-count ORFs are filtered before fitting). Per-contrast results land under `<outdir>/dte/orf_level/<method>/`; the CDS-aggregated gene matrix under `<outdir>/dte/gene_level_cds_aggregated/`. CLI overrides for the ORF-level deltaTE and DOTSeq fits are exposed via `--extra_orf_dte_args` and `--extra_dotseq_args`. Row-independence caveat for Tier 2 (ORFs sharing a gene-level RNA-seq denominator are perfectly correlated after the join) is documented in `docs/usage.md` ([@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [#133](https://github.com/nf-core/riboseq/pull/133) - Improve anota2seq documentation: fix typos in contrast file docs, document available `extra_anota2seq_run_args` options ([@pinin4fjords](https://github.com/pinin4fjords))
- [#145](https://github.com/nf-core/riboseq/pull/145) - Filter to primary alignments before UMI-tools dedup ([@pinin4fjords](https://github.com/pinin4fjords))
- [#147](https://github.com/nf-core/riboseq/pull/147) - Fix RiboCode error handling and P-site detection failures ([@pinin4fjords](https://github.com/pinin4fjords))
- [#148](https://github.com/nf-core/riboseq/pull/148) - Update all nf-core modules/subworkflows to latest versions, migrate to topic-based version reporting, fix publishing paths for new subworkflow nesting ([@pinin4fjords](https://github.com/pinin4fjords))
- [#155](https://github.com/nf-core/riboseq/pull/155) - Drop `if (params.skip_X)` wrappers in `conf/modules.config` so the Nextflow 26.04 v2 strict config parser accepts the file ([@pinin4fjords](https://github.com/pinin4fjords))

### `Changed`

- [#168](https://github.com/nf-core/riboseq/issues/168) - Extend `--translational_efficiency_method` to cover the Tier 2 ORF-level DTE path. anota2seq, deltaTE and DOTSeq all run at ORF resolution against a common method-agnostic combined count matrix produced from `orf_psite_counts.tsv` joined with the gene-level Salmon RNA-seq matrix via `orf_to_gene.tsv`. DOTSeq additionally fits its package-specific Differential ORF Usage (DOU) contrast - a per-gene beta-binomial GLM modelling whether each ORF gains or loses share of its parent gene's ribosome occupancy across conditions (no anota2seq / deltate counterpart). Selecting `dotseq` requires `--extended_orf_analysis true` and at least one ORF caller. ORF-level results land under `<outdir>/dte/orf_level/<method>/` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#129](https://github.com/nf-core/riboseq/pull/129) - Bump pipeline version to 1.3.0dev ([@iraiosub](https://github.com/iraiosub))
- [#152](https://github.com/nf-core/riboseq/pull/152) - Update ribodetector module to 0.3.3 (nf-core/modules#11131) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Update nf-core/multiqc module to 1.34 ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Align FastQC/SortMeRNA general-stats handling with `nf-core/rnaseq` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#160](https://github.com/nf-core/riboseq/issues/160) - **Breaking default change:** flip the default `--te_quantification_method` from `alignment` (STAR + Salmon) to `plastid_psite` (in-frame P-site counts from plastid). Salmon's coverage-uniformity, fragment-length and multi-mapping assumptions are inappropriate for short, length-constrained Ribo-seq footprints; in-frame P-site counts are the scientifically correct quantity. The two methods produce fundamentally different per-gene count values (Salmon TPM-derived pseudo-counts vs raw integer in-frame P-site sums), so re-running an existing cohort on a newer pipeline version with the default will not reproduce the previous count matrix. Users who need backward-compatible Salmon-style counts must now set `--te_quantification_method alignment` explicitly. ([@pinin4fjords](https://github.com/pinin4fjords))
- [#163](https://github.com/nf-core/riboseq/issues/163) - Demote Ribotricer from the default ORF caller set to opt-in: replace `--skip_ribotricer` (default `false`) with `--run_ribotricer` (default `false`), warn at runtime when enabled, and exclude its score column from cross-caller rank aggregation. The default caller set is now Ribo-TISH + RiboCode; agreement logic is parameterised on the runtime-enabled caller set ([@pinin4fjords](https://github.com/pinin4fjords))

### `Parameters`

| Old parameter       | New parameter                     |
| ------------------- | --------------------------------- |
| `--skip_ribotricer` | `--run_ribotricer`                |
|                     | `--canonical_gtf`                 |
|                     | `--skip_stringtie`                |
|                     | `--extra_stringtie_args`          |
|                     | `--extra_stringtie_merge_args`    |
|                     | `--stringtie_class_codes`         |
|                     | `--stringtie_ribo_fallback_args`  |
|                     | `--novel_gtf`                     |
|                     | `--rrna_blacklist`                |
|                     | `--extended_orf_analysis`         |
|                     | `--run_rpbp`                      |
|                     | `--extra_rpbp_preparegenome_args` |
|                     | `--extra_rpbp_predictorfs_args`   |
|                     | `--run_price`                     |
|                     | `--extra_price_indexgenome_args`  |
|                     | `--extra_price_price_args`        |
|                     | `--extra_orf_dte_args`            |
|                     | `--extra_dotseq_args`             |

### `Dependencies`

| Dependency         | Old version | New version |
| ------------------ | ----------- | ----------- |
| `MultiQC`          | 1.32        | 1.33        |
| `plastid`          |             | 0.6.1       |
| `bedtools`         |             | 2.31.1      |
| `bedGraphToBigWig` |             | 469         |
| `AGAT`             |             | 1.6.1       |
| `StringTie`        |             | 2.2.3       |

## v1.2.0 - 2025-12-03

### `Added`

- [#113](https://github.com/nf-core/riboseq/pull/113) - Add P-site identification with riboWaltz ([@iraiosub](https://github.com/iraiosub), reviewed by [@JackCurragh](https://github.com/JackCurragh))

### `Parameters`

| Old parameter | New parameter            |
| ------------- | ------------------------ |
|               | `--skip_ribowaltz`       |
|               | `--extra_ribowaltz_args` |
|               | `--fastp_merge`          |

### `Changed`

- [#103](https://github.com/nf-core/riboseq/pull/103) - Updated the JSON schema to make input validation stricter([@andreirie](https://github.com/andreirie))
- [#118](https://github.com/nf-core/riboseq/pull/118) - Template update for nf-core/tools v3.5.1 ([@nf-corebot](https://github.com/nf-corebot))
- [#111](https://github.com/nf-core/riboseq/pull/111) - Template update for nf-core/tools v3.4.1 ([@nf-corebot](https://github.com/nf-corebot), [@maxulysse](https://github.com/maxulysse), reviewed by [@mashehu](https://github.com/mashehu))
- [#115](https://github.com/nf-core/riboseq/pull/115) - Prerelease changes v1.2.0 ([@iraiosub](https://github.com/iraiosub), reviewed by [@JackCurragh](https://github.com/JackCurragh))
- [#117](https://github.com/nf-core/riboseq/pull/117) - Update modules: `fq/lint`, `fq/subsample` and `sortmerna`([@iraiosub](https://github.com/iraiosub), reviewed by [@mashehu](https://github.com/mashehu))
- [#120](https://github.com/nf-core/riboseq/pull/120) - Bump Nextflow minimum version to 25.04.8 ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@iraiosub](https://github.com/iraiosub))
- [#121](https://github.com/nf-core/riboseq/pull/121) - Remove conda from CI due to upstream ribotish Python 3.14 incompatibility ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@iraiosub](https://github.com/iraiosub))

### `Fixed`

- [#114](https://github.com/nf-core/riboseq/pull/114) - Fix order of steps in metro map ([#109](https://github.com/nf-core/riboseq/issues/109)): UMI-tools extract is now correctly placed before adaptor trimming ([@iraiosub](https://github.com/iraiosub), reviewed by [@JackCurragh](https://github.com/JackCurragh))
- [#117](https://github.com/nf-core/riboseq/pull/117) - Update `ribowaltz` module containing a conda dependencies fix and add missing versions to MultiQC ([@iraiosub](https://github.com/iraiosub), reviewed by [@mashehu](https://github.com/mashehu))
- [#124](https://github.com/nf-core/riboseq/pull/124) - Placed `ribowaltz` versions inside conditional ([@iraiosub](https://github.com/iraiosub), reviewed by [@pinin4fjords](https://github.com/pinin4fjords))
- [#126](https://github.com/nf-core/riboseq/pull/126) - Update `fastq_qc_trim_filter_setstrandedness` subworkflow to ensure compatibility with the `fastq_fastqc_umitools_fastp` subworkflow ([@iraiosub](https://github.com/iraiosub))

### `Dependencies`

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `fastp`     | 0.23.4      | 1.0.1       |
| `MultiQC`   | 1.27        | 1.32        |
| `umitools`  | 1.1.5       | 1.1.6       |
| `SortMeRNA` | 4.3.6       | 4.3.7       |
| `Nextflow`  | 24.04.2     | 25.04.8     |

### `Deprecated`

## v1.1.0 - 2025-01-30

### `Added`

### `Changed`

- [#61](https://github.com/nf-core/riboseq/pull/61) - Update Metro Map ([@maxulysse](https://github.com/maxulysse), reviewed by [@drpatelh](https://github.com/drpatelh))
- [#71](https://github.com/nf-core/riboseq/pull/71) - Template update for nf-core/tools v3.0.2 ([@nf-corebot](https://github.com/nf-corebot), ([@maxulysse](https://github.com/maxulysse), reviewed by [@JackCurragh](https://github.com/JackCurragh)), [@FelixKrueger](https://github.com/FelixKrueger))
- [#73](https://github.com/nf-core/riboseq/pull/73) - Pipeline level snapshots with nf-test (([@maxulysse](https://github.com/maxulysse), reviewed by [@pinin4fjords](https://github.com/pinin4fjords)))
- [#77](https://github.com/nf-core/riboseq/pull/77) - Update `RIBOTRICER_PREPAREORFS` to increase resource allocation ([@iraiosub](https://github.com/iraiosub))
- [#83](https://github.com/nf-core/riboseq/pull/83) - Fix skip_ribotish ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#85](https://github.com/nf-core/riboseq/pull/85) - Rationalise sorted bam/bai publishing ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#86](https://github.com/nf-core/riboseq/pull/86) - Important! Template update for nf-core/tools v3.2.0 ([@nf-core-bot](https://github.com/nf-core-bot), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [#87](https://github.com/nf-core/riboseq/pull/87) - Bump versions pre-release to 1.1.0 ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#92](https://github.com/nf-core/riboseq/pull/92) - Bump anota2seq for ordering fix ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#93](https://github.com/nf-core/riboseq/pull/93) - Bump anota2seq for dollar fix ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#94](https://github.com/nf-core/riboseq/pull/94) - Fix value channel for multi-contrast case ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#96](https://github.com/nf-core/riboseq/pull/96) - Fix minor linting issue for release ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#97](https://github.com/nf-core/riboseq/pull/97) - Remove the on_pull_request_target trigger from the download test ([@pinin4fjords](https://github.com/pinin4fjords))
- [#98](https://github.com/nf-core/riboseq/pull/98) - Bump gunzip due to release CI failure ([@pinin4fjords](https://github.com/pinin4fjords))
- [#99](https://github.com/nf-core/riboseq/pull/99) - Remove conda from release CI ([@pinin4fjords](https://github.com/pinin4fjords))
- [#100](https://github.com/nf-core/riboseq/pull/100) - Fix gunzip in snapshot ([@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [#60](https://github.com/nf-core/riboseq/pull/60) - Pass empty value to samtools sort in UMI handling branch ([@JackCurragh](https://github.com/JackCurragh), reviewed by [@FelixKrueger](https://github.com/FelixKrueger), [@maxulysse](https://github.com/maxulysse))
- [#61](https://github.com/nf-core/riboseq/pull/61) - Update subworkflow `utils_nfcore_pipeline` ([@maxulysse](https://github.com/maxulysse), reviewed by [@drpatelh](https://github.com/drpatelh))
- [#75](https://github.com/nf-core/riboseq/pull/75) - UMI fixes: solve deduplication issue and update input handling for Salmon ([@iraiosub](https://github.com/iraiosub), reviewed by [@FelixKrueger](https://github.com/FelixKrueger), [@pinin4fjords](https://github.com/pinin4fjords))
- [#79](https://github.com/nf-core/riboseq/pull/75) - Move UMI handling to subworkflow, update modules and subworkflows, deal with docs and config fallout ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))
- [#90](https://github.com/nf-core/riboseq/pull/90) - --subset_to_contrast_samples must be true for anota2seq ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `MultiQC`  | 1.21        | 1.25.1      |

### `Parameters`

| Old parameter                        | New parameter         |
| ------------------------------------ | --------------------- |
|                                      | `--help_full`         |
|                                      | `--show_hidden`       |
|                                      | `--skip_linting`      |
|                                      | `--extra_fqlint_args` |
|                                      | `--umi_dedup_tool`    |
| `--validationFailUnrecognisedParams` |                       |
| `--validationLenientMode`            |                       |
| `--validationSchemaIgnoreParams`     |                       |
| `--validationShowHiddenParams`       |                       |

### `Deprecated`

## v1.0.1 - 2024-04-17

### `Added`

- [#53](https://github.com/nf-core/riboseq/pull/53) - Bump to v1.0.1 and set Zenodo ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#54](https://github.com/nf-core/riboseq/pull/54) - Add legal logos and first metro map ([@FelixKrueger](https://github.com/FelixKrueger), review by [@maxulysse](https://github.com/maxulysse))

### `Fixed`

- [#57](https://github.com/nf-core/riboseq/pull/57) - Minor spacing changes to logo svg ([@JackCurragh](https://github.com/JackCurragh), reviewed by [@FelixKrueger](https://github.com/FelixKrueger))

## v1.0.0 - 2024-04-12

Initial release of nf-core/riboseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#3](https://github.com/nf-core/riboseq/pull/3) - Re-initialise base template ([@maxulysse](https://github.com/maxulysse), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [#4](https://github.com/nf-core/riboseq/pull/4) - Initialise testing an base template functionality ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#8](https://github.com/nf-core/riboseq/pull/8) - Preprocessing from rnaseq ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse) review by [@adamrtalbot](https://github.com/adamrtalbot))
- [#10](https://github.com/nf-core/riboseq/pull/10) - Take preprocessing from nf-core ([@pinin4fjords](https://github.com/pinin4fjords), review by [@adamrtalbot](https://github.com/adamrtalbot))
- [#12](https://github.com/nf-core/riboseq/pull/12) - Add alignment via STAR + postprocessing (([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#35](https://github.com/nf-core/riboseq/pull/35) - Sortmerna: index once ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#40](https://github.com/nf-core/riboseq/pull/40) - Ribotricer orf prediction ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#42](https://github.com/nf-core/riboseq/pull/42) - Add alignment based quantification with Salmon ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#43](https://github.com/nf-core/riboseq/pull/43) - Add translational efficiency analysis with anota2seq ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))

### `Changed`

- [#9](https://github.com/nf-core/riboseq/pull/9) - Important! Template update for nf-core/tools v2.12 ([nf-core-bot](https://github.com/nf-core-bot), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [#32](https://github.com/nf-core/riboseq/pull/32) - Nf core template merge 2.13 (manual) ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse), [@adamrtalbot](https://github.com/adamrtalbot))
- [#38](https://github.com/nf-core/riboseq/pull/38) - Important! Template update for nf-core/tools v2.13.1 ([@nf-core-bot](https://github.com/nf-core-bot), [@pinin4fjords](https://github.com/pinin4fjords))
- [#46](https://github.com/nf-core/riboseq/pull/46) - Prerelease changes v1.0.0 ([@pinin4fjords](https://github.com/pinin4fjords), review by [@FelixKrueger](https://github.com/FelixKrueger))
- [#51](https://github.com/nf-core/riboseq/pull/51) - Change to custom logo ([@JackCurragh](https://github.com/jackcurragh), review by [@FelixKrueger](https://github.com/FelixKrueger))

### `Fixed`

- [#5](https://github.com/nf-core/riboseq/pull/5) - Fix linting ([@maxulysse](https://github.com/maxulysse), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [#34](https://github.com/nf-core/riboseq/pull/34) - Fix order of preprocessing steps ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#36](https://github.com/nf-core/riboseq/pull/36) - Bump bbsplit module to prevent index overwrites ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#44](https://github.com/nf-core/riboseq/pull/44) - Fix lack of fastqc in multiqc report ([@pinin4fjords](https://github.com/pinin4fjords), review by [@mashehu](https://github.com/mashehu))
- [#45](https://github.com/nf-core/riboseq/pull/45) - Update CI from rnaseq, strip unused rnaseq components ([@pinin4fjords](https://github.com/pinin4fjords), review by [@jfy133](https://github.com/jfy133))
- [#48](https://github.com/nf-core/riboseq/pull/48) - Remove stub option from download in CI ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [#49](https://github.com/nf-core/riboseq/pull/49) - Fix CI ([@pinin4fjords](https://github.com/pinin4fjords), review by [@adamrtalbot](https://github.com/adamrtalbot))
- [#50](https://github.com/nf-core/riboseq/pull/50) - V1.0.0 release review fixes ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))

### `Dependencies`

### `Deprecated`
