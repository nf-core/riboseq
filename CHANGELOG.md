# nf-core/riboseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev

### `Added`

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
- [#161](https://github.com/nf-core/riboseq/issues/161) - Add a one-transcript-per-gene canonical annotation backbone via the new `--canonical_gtf` parameter, used for ORF calling, riboWaltz P-site calibration, plastid P-site quantification and the translational-efficiency analysis; falls back to AGAT longest-isoform extraction from `--gtf` when not supplied. The full `--gtf` is retained for genome-guided alignment ([@pinin4fjords](https://github.com/pinin4fjords))
- [#179](https://github.com/nf-core/riboseq/issues/179) - Wire the optional secondary reference annotation (`-a`) through to ribotish/predict via a dedicated 4-tuple `ch_fasta_gtf_for_ribotish` channel ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Novel transcript discovery: StringTie reference-guided assembly (RNA-seq BAMs preferred, Ribo-seq fallback with tightened defaults) or user-supplied GTF via `--novel_gtf`. Output passes through gffcompare class-code filter (default `u` intergenic only), optional strand-aware rRNA/repeat blacklist (`--rrna_blacklist`), and is concatenated with the canonical backbone into `<outdir>/stringtie/hybrid_reference.gtf`. The hybrid GTF is exposed on the `hybrid_gtf` workflow channel (equals canonical when no novel source is configured) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#165](https://github.com/nf-core/riboseq/issues/165) - Add `--extended_orf_analysis` (default `false`) which routes the hybrid GTF into the genome-BAM ORF callers (Ribo-TISH `predict`, Ribotricer `prepare-orfs`) so novel intergenic ORFs are within scope. RiboCode and transcriptome-BAM consumers stay on the canonical backbone. When the flag is set without a novel-transcript source, the pipeline warns and falls back to canonical ([@pinin4fjords](https://github.com/pinin4fjords))
- [#171](https://github.com/nf-core/riboseq/issues/171) - Second STAR pass against a hybrid transcriptome (gffread-extracted from `hybrid_reference.gtf`) when `--extended_orf_analysis true`. RiboCode is rewired to consume the resulting hybrid transcriptome BAM and the hybrid GTF, bringing novel intergenic transcripts into scope. Ribo-seq samples only - RNA-seq and TI-seq are not consumed by RiboCode. riboWaltz, plastid and Salmon-based quantification continue on the canonical reference BAM by design. Hybrid outputs land under `<outdir>/hybrid_star/` ([@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [#133](https://github.com/nf-core/riboseq/pull/133) - Improve anota2seq documentation: fix typos in contrast file docs, document available `extra_anota2seq_run_args` options ([@pinin4fjords](https://github.com/pinin4fjords))
- [#145](https://github.com/nf-core/riboseq/pull/145) - Filter to primary alignments before UMI-tools dedup ([@pinin4fjords](https://github.com/pinin4fjords))
- [#147](https://github.com/nf-core/riboseq/pull/147) - Fix RiboCode error handling and P-site detection failures ([@pinin4fjords](https://github.com/pinin4fjords))
- [#148](https://github.com/nf-core/riboseq/pull/148) - Update all nf-core modules/subworkflows to latest versions, migrate to topic-based version reporting, fix publishing paths for new subworkflow nesting ([@pinin4fjords](https://github.com/pinin4fjords))
- [#155](https://github.com/nf-core/riboseq/pull/155) - Drop `if (params.skip_X)` wrappers in `conf/modules.config` so the Nextflow 26.04 v2 strict config parser accepts the file ([@pinin4fjords](https://github.com/pinin4fjords))

### `Changed`

- [#129](https://github.com/nf-core/riboseq/pull/129) - Bump pipeline version to 1.3.0dev ([@iraiosub](https://github.com/iraiosub))
- [#152](https://github.com/nf-core/riboseq/pull/152) - Update ribodetector module to 0.3.3 (nf-core/modules#11131) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Update nf-core/multiqc module to 1.34 ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Align FastQC/SortMeRNA general-stats handling with `nf-core/rnaseq` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#160](https://github.com/nf-core/riboseq/issues/160) - **Breaking default change:** flip the default `--te_quantification_method` from `alignment` (STAR + Salmon) to `plastid_psite` (in-frame P-site counts from plastid). Salmon's coverage-uniformity, fragment-length and multi-mapping assumptions are inappropriate for short, length-constrained Ribo-seq footprints; in-frame P-site counts are the scientifically correct quantity. The two methods produce fundamentally different per-gene count values (Salmon TPM-derived pseudo-counts vs raw integer in-frame P-site sums), so re-running an existing cohort on a newer pipeline version with the default will not reproduce the previous count matrix. Users who need backward-compatible Salmon-style counts must now set `--te_quantification_method alignment` explicitly. ([@pinin4fjords](https://github.com/pinin4fjords))
- [#162](https://github.com/nf-core/riboseq/issues/162) - Document that Ribo-TISH `quality` continues to consume the canonical annotation (P-site offset calibration uses CDS-bearing canonical transcripts; mixing CDS-absent novel transcripts would degrade calibration without adding diagnostic signal) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#163](https://github.com/nf-core/riboseq/issues/163) - Demote Ribotricer from the default ORF caller set to opt-in: replace `--skip_ribotricer` (default `false`) with `--run_ribotricer` (default `false`), warn at runtime when enabled, and exclude its score column from cross-caller rank aggregation. The default caller set is now Ribo-TISH + RiboCode; agreement logic is parameterised on the runtime-enabled caller set ([@pinin4fjords](https://github.com/pinin4fjords))

### `Parameters`

| Old parameter       | New parameter                            |
| ------------------- | ---------------------------------------- |
| `--skip_ribotricer` | `--run_ribotricer`                       |
|                     | `--ribo_removal_tool`                    |
|                     | `--skip_ribocode`                        |
|                     | `--extra_ribocode_gtfupdate_args`        |
|                     | `--extra_ribocode_prepare_args`          |
|                     | `--extra_ribocode_metaplots_args`        |
|                     | `--extra_ribocode_ribocode_args`         |
|                     | `--translational_efficiency_method`      |
|                     | `--extra_deltate_args`                   |
|                     | `--te_lfc_threshold`                     |
|                     | `--rna_lfc_threshold`                    |
|                     | `--ribo_lfc_threshold`                   |
|                     | `--equalise_read_lengths`                |
|                     | `--equalise_read_lengths_target`         |
|                     | `--skip_plastid`                         |
|                     | `--plastid_min_length`                   |
|                     | `--plastid_max_length`                   |
|                     | `--plastid_default_psite_offset`         |
|                     | `--extra_plastid_metagene_generate_args` |
|                     | `--extra_plastid_psite_args`             |
|                     | `--extra_plastid_make_wiggle_args`       |
|                     | `--skip_coverage_tracks`                 |
|                     | `--ribodetector_chunk_size`              |
|                     | `--canonical_gtf`                        |
|                     | `--skip_stringtie`                       |
|                     | `--novel_gtf`                            |
|                     | `--stringtie_class_codes`                |
|                     | `--rrna_blacklist`                       |
|                     | `--extra_stringtie_args`                 |
|                     | `--extra_stringtie_merge_args`           |
|                     | `--stringtie_ribo_fallback_args`         |
|                     | `--extended_orf_analysis`                |

### `Dependencies`

| Dependency         | Old version | New version |
| ------------------ | ----------- | ----------- |
| `MultiQC`          | 1.32        | 1.33        |
| `plastid`          |             | 0.6.1       |
| `bedtools`         |             | 2.31.1      |
| `bedGraphToBigWig` |             | 469         |
| `AGAT`             |             | 1.6.1       |

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
