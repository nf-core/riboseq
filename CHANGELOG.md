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
- [#154](https://github.com/nf-core/riboseq/pull/154) - Add riboWaltz QC plots to the MultiQC report via a new MultiQC riboWaltz module ([@pinin4fjords](https://github.com/pinin4fjords))
- [#161](https://github.com/nf-core/riboseq/issues/161) - Add a one-transcript-per-gene canonical annotation backbone via `--canonical_gtf` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#179](https://github.com/nf-core/riboseq/issues/179) - Wire an optional secondary reference annotation (`-a`) through to ribotish/predict ([@pinin4fjords](https://github.com/pinin4fjords))
- [#164](https://github.com/nf-core/riboseq/issues/164) - Add novel transcript discovery via StringTie or `--novel_gtf`, merged with the canonical backbone into a hybrid reference GTF ([@pinin4fjords](https://github.com/pinin4fjords))
- [#165](https://github.com/nf-core/riboseq/issues/165) - Add `--extended_orf_analysis` to route the hybrid GTF into the genome-BAM ORF callers ([@pinin4fjords](https://github.com/pinin4fjords))
- [#171](https://github.com/nf-core/riboseq/issues/171) - Add a second STAR pass against the hybrid transcriptome so RiboCode can call novel-transcript ORFs under `--extended_orf_analysis` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#170](https://github.com/nf-core/riboseq/issues/170) - Add PRICE as an opt-in ORF caller via `--run_price` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#169](https://github.com/nf-core/riboseq/issues/169) - Add Rp-Bp as an opt-in ORF caller via `--run_rpbp` ([@pinin4fjords](https://github.com/pinin4fjords))
- [#167](https://github.com/nf-core/riboseq/issues/167) - Add a cross-caller cohort ORF catalogue under `--extended_orf_analysis`, with smORF peptide collapse (`--skip_orf_collapse`) and a consensus view (`--orf_min_callers`/`--orf_min_samples`) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#166](https://github.com/nf-core/riboseq/issues/166) - Add per-ORF in-frame P-site quantification, emitting an ORF x sample count matrix ([@pinin4fjords](https://github.com/pinin4fjords))
- [#168](https://github.com/nf-core/riboseq/issues/168) - Add ORF-level differential translation analysis (anota2seq / deltaTE / DOTSeq) on top of the gene-level DTE ([@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [#133](https://github.com/nf-core/riboseq/pull/133) - Improve anota2seq documentation: fix typos in contrast file docs, document available `extra_anota2seq_run_args` options ([@pinin4fjords](https://github.com/pinin4fjords))
- [#145](https://github.com/nf-core/riboseq/pull/145) - Filter to primary alignments before UMI-tools dedup ([@pinin4fjords](https://github.com/pinin4fjords))
- [#147](https://github.com/nf-core/riboseq/pull/147) - Fix RiboCode error handling and P-site detection failures ([@pinin4fjords](https://github.com/pinin4fjords))
- [#148](https://github.com/nf-core/riboseq/pull/148) - Update all nf-core modules/subworkflows to latest versions, migrate to topic-based version reporting, fix publishing paths for new subworkflow nesting ([@pinin4fjords](https://github.com/pinin4fjords))
- [#155](https://github.com/nf-core/riboseq/pull/155) - Drop `if (params.skip_X)` wrappers in `conf/modules.config` so the Nextflow 26.04 v2 strict config parser accepts the file ([@pinin4fjords](https://github.com/pinin4fjords))
- [#193](https://github.com/nf-core/riboseq/pull/193) - Re-pin `custom/orfnormalise` (nf-core/modules#12173) to fix ribotish `GenomePos` coordinate parsing (0-based, was treated as 1-based), which shifted ribotish-derived catalogue ORFs out of frame ([@pinin4fjords](https://github.com/pinin4fjords))
- [#194](https://github.com/nf-core/riboseq/pull/194) - Update `dotseq/dotseq` (nf-core/modules#12230) so the DOTSeq heatmap actually renders, and pass `--orf_type_main_value canonical_cds` to pair the heatmap against riboseq's main-ORF label ([@pinin4fjords](https://github.com/pinin4fjords))
- [#204](https://github.com/nf-core/riboseq/pull/204) - Include ORFs on novel transcripts in ORF-level DTE, by quantifying the RNA-seq denominator against the full reference transcriptome augmented with the novel intergenic transcripts so novel genes gain a host-gene row ([@pinin4fjords](https://github.com/pinin4fjords))

### `Changed`

- [#129](https://github.com/nf-core/riboseq/pull/129) - Bump pipeline version to 1.3.0dev ([@iraiosub](https://github.com/iraiosub))
- [#152](https://github.com/nf-core/riboseq/pull/152) - Update ribodetector module to 0.3.3 (nf-core/modules#11131) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Add riboWaltz QC plots to the MultiQC report via a new MultiQC riboWaltz module ([@pinin4fjords](https://github.com/pinin4fjords))
- [#154](https://github.com/nf-core/riboseq/pull/154) - Add riboWaltz QC plots to the MultiQC report via a new MultiQC riboWaltz module ([@pinin4fjords](https://github.com/pinin4fjords))
- [#160](https://github.com/nf-core/riboseq/issues/160) - Change the default `--te_quantification_method` from `alignment` (STAR + Salmon) to `plastid_psite` (in-frame P-site counts). **Breaking:** per-gene counts are not comparable to Salmon-based runs; pass `--te_quantification_method alignment` to restore the previous behaviour. ([@pinin4fjords](https://github.com/pinin4fjords))
- [#162](https://github.com/nf-core/riboseq/issues/162) - Document that Ribo-TISH `quality` continues to consume the canonical annotation (P-site offset calibration uses CDS-bearing canonical transcripts; mixing CDS-absent novel transcripts would degrade calibration without adding diagnostic signal) ([@pinin4fjords](https://github.com/pinin4fjords))
- [#163](https://github.com/nf-core/riboseq/issues/163) - Demote Ribotricer from the default ORF caller set to opt-in: replace `--skip_ribotricer` (default `false`) with `--run_ribotricer` (default `false`), warn at runtime when enabled, and exclude its score column from cross-caller rank aggregation. The default caller set is now Ribo-TISH + RiboCode; agreement logic is parameterised on the runtime-enabled caller set ([@pinin4fjords](https://github.com/pinin4fjords))
- [#194](https://github.com/nf-core/riboseq/pull/194) - Update all nf-core modules and subworkflows to their latest revisions and reconcile the pipeline with the changed component interfaces ([@pinin4fjords](https://github.com/pinin4fjords))
- [#194](https://github.com/nf-core/riboseq/pull/194) - Raise the minimum Nextflow version to `25.10.0`, required by nf-schema 2.7.2 ([@pinin4fjords](https://github.com/pinin4fjords))
- [#156](https://github.com/nf-core/riboseq/pull/156) - Template update for nf-core/tools v4.0.2 ([@nf-core-bot](https://github.com/nf-core-bot), [@pinin4fjords](https://github.com/pinin4fjords))

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
|                     | `--gffcompare_class_codes`               |
|                     | `--rrna_blacklist`                       |
|                     | `--extra_stringtie_args`                 |
|                     | `--extra_stringtie_merge_args`           |
|                     | `--extended_orf_analysis`                |
|                     | `--run_price`                            |
|                     | `--extra_price_indexgenome_args`         |
|                     | `--extra_price_price_args`               |
|                     | `--run_rpbp`                             |
|                     | `--extra_rpbp_preparegenome_args`        |
|                     | `--extra_rpbp_predictorfs_args`          |
|                     | `--extra_orf_deltate_args`               |
|                     | `--extra_orf_anota2seq_run_args`         |
|                     | `--extra_dotseq_args`                    |

### `Dependencies`

| Dependency         | Old version | New version |
| ------------------ | ----------- | ----------- |
| `MultiQC`          | 1.32        | 1.35        |
| `nf-schema`        | 2.5.1       | 2.7.2       |
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
