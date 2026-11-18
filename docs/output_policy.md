# Output policy

This document defines what this fork of nf-core/riboseq saves, what it exports for analysis, and what is kept only for troubleshooting. For a description of each output file, see [output.md](output.md).

> Status: implemented in the pipeline; not yet in a tagged release.

## Generated vs exported

The pipeline keeps nf-core's `results/` layout and publishes it on CAVATICA as one directory.

- **What gets saved at all** is controlled by `publishDir` in `conf/modules.config`. Anything not needed is not published, rather than published and deleted later.
- **What gets delivered** is the files listed in `export_manifest.txt`, plus the whole `pipeline_info/` folder.
- **Delivery is manual for now.** Whoever delivers a run copies those files from the CAVATICA task outputs to the delivery location, using the manifest as the checklist. No pipeline step copies or deletes files.
- Files not delivered can be deleted from the CAVATICA project after QC sign-off (MultiQC reviewed). Re-run the pipeline if they are needed later.

The nf-core layout is kept, rather than flattened, so that upstream updates can still be merged cleanly.

### How the manifest is written

When a run succeeds, the pipeline's completion handler (`writeExportManifest` in `subworkflows/local/utils_nfcore_riboseq_pipeline/main.nf`) writes `<outdir>/export_manifest.txt`. It lists every published file that matches `params.export_paths`, one path per line relative to `outdir`, sorted.

- Patterns are globs relative to `outdir`. `*` stays within one folder; `**` crosses folders. A pattern starting with `!` removes matching files.
- A pattern that matches nothing is logged as a warning, not an error, because some outputs exist only for some runs (for example, translational efficiency needs contrasts).
- A failed run writes no manifest.
- `pipeline_info/` is delivered as a whole folder, because Nextflow writes its HTML reports after the manifest.

## Alignments

All saved alignments are CRAM (`alignment_output_format = 'cram'`). BAMs exist only temporarily inside a run.

- Genome and transcriptome alignments are both saved as CRAM, so P-site and quantification steps can be re-run without realigning.
- A CRAM can only be read with the reference it was made with, so `genome/cram_reference/genome.*` and `transcriptome.*` are exported with every run.

## Not saved

| Output | Why |
|--------|-----|
| Intermediate alignment BAMs | `save_align_intermeds = false`. In CRAM mode only CRAMs are published anyway. |
| `genome/` (indexes, sortmerna, FASTA copies), except `genome/cram_reference/` | Rebuilt from the reference (`save_reference = false`). |
| `orf_predictions/ribocode/annotation/` | Rebuilt from the reference; identical across runs that use the same genome. |

## Default export list

`params.export_paths`, relative to `outdir`:

| Path | Why |
|------|-----|
| `alignment/star/sorted/*.sorted.cram` (+ `.crai`) | Genome and transcriptome alignments |
| `genome/cram_reference/genome.*`, `transcriptome.*` | References needed to read the CRAMs |
| `multiqc/**` | QC report |
| `coverage/*.bigWig` | Tracks for IGV or the UCSC genome browser |
| `psites/ribowaltz/ribowaltz_qc/*` | P-site QC plots and summary tables (not the per-read P-site TSV.gz files) |
| `quantification/salmon/salmon.merged.*` | Gene and transcript counts, TPM and lengths (TSV), plus the SummarizedExperiment RDS files |
| `quantification/inframe_psite/gene_counts.tsv` | In-frame P-site gene counts |
| `quantification/orf_level/orf_psite_counts.tsv` | ORF count matrix |
| `orf_predictions/catalogue/**` | ORF catalogue: consensus, summary (`cohort_orfs.tsv`, per-sample ORF tables) |
| `orf_predictions/ribotish_all/allsamples_*` | Joint Ribo-TISH predictions |
| `translational_efficiency/**.tsv`, `**.rds` | deltaTE, anota2seq and DOTSeq results with gene symbols |
| `transcript_assembly/stringtie/*.gtf`, excluding `!…/*.denovo.transcripts.gtf` | Merged transcript models that the ORF catalogue refers to |

Plus `pipeline_info/` (parameters, software versions, execution trace and reports), delivered as a folder.

RDS files are exported only for genuine R objects (SummarizedExperiment, DESeqDataSet, Anota2seq, DOTSeq). Each has TSV tables beside it for readers who don't use R.

## Generated but not exported (troubleshooting only)

- `preprocessing/` (trimming, rRNA removal, FastQC, read-length statistics). These are summarized in MultiQC.
- `alignment/*/log/`, `samtools_stats/`, `picard_metrics/`
- Hybrid-STAR alignments and their reference (`alignment/hybrid_star/`, `genome/cram_reference/hybrid_transcriptome.*`), produced only with `extended_orf_analysis`
- Per-tool ORF predictions (`orf_predictions/ribotish/`, `ribocode/`, `ribotricer/`, `price/`, `rpbp/`)
- The per-read riboWaltz P-site tables, `psites/plastid/`, and `riboseq_qc/`
- The per-sample StringTie GTFs and `gffcompare/`
- Translational-efficiency plots

## Merging across samples

Merging happens inside a single run, over the samplesheet: the Salmon merge, the ORF count matrix, the ORF consensus, and the translational-efficiency contrasts. There is no separate merge workflow. When new samples arrive, the whole cohort is re-run.

### Alternative for discussion: separate per-sample and cohort stages

The pipeline could instead be split into two workflows:

1. A **per-sample** workflow (trimming, rRNA removal, alignment, P-site calling, per-sample quantification), run once per sample as it arrives.
2. A **cohort** workflow that takes the per-sample outputs and runs the Salmon merge, joint ORF calling, the ORF consensus and count matrix, and the translational-efficiency contrasts.

**For:** new samples don't require realigning old ones, and cohort-level analysis can be re-run cheaply with different contrasts.

**Against:**
- It diverges substantially from upstream nf-core/riboseq, making merges harder.
- Joint ORF calling (Ribo-TISH `allsamples_*`) and the ORF consensus depend on all samples, so the cohort stage remains the expensive part.
- It needs a new CAVATICA app and a merge ticket.

## Versioning

Releases use semantic versioning (`vMAJOR.MINOR.PATCH`), consistent with other childrens-bti repositories. The first release is `v1.0.0`, based on nf-core/riboseq 2.0.0. `manifest.version` in `nextflow.config` stays at the upstream version so nf-core tooling and upstream merges keep working. The existing `1.0.0`–`1.2.0` tags (no `v`) are upstream nf-core releases. Containers are pinned per module, as is standard for nf-core.
