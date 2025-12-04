# nf-core/riboseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev

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
