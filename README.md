<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-riboseq_logo_dark.png">
    <img alt="nf-core/riboseq" src="docs/images/nf-core-riboseq_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/riboseq)
[![GitHub Actions CI Status](https://github.com/nf-core/riboseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/riboseq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/riboseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/riboseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/riboseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.10966364-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.10966364)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.8-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/riboseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23riboseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/riboseq)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/riboseq** is a bioinformatics pipeline for analysis of Ribo-seq data. It borrows heavily from nf-core/rnaseq in the preprocessing stages:

![nf-core/riboseq metro map](docs/images/nf-core-riboseq_metro_map.png)

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Sub-sample FastQ files and auto-infer strandedness ([`fq`](https://github.com/stjude-rust-labs/fq), [`Salmon`](https://combine-lab.github.io/salmon/))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
5. Adapter and quality trimming ([`Trim Galore!`](https://github.com/FelixKrueger/TrimGalore))
6. Removal of genome contaminants ([`BBSplit`](http://seqanswers.com/forums/showthread.php?t=41288))
7. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
8. Genome alignment of reads, outputting both genome and transcriptome alignments with [`STAR`](https://github.com/alexdobin/STAR)
9. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
10. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))

Differences occur in the downstream analysis steps. Currently these specialist steps are:

1. (default) Reduce the reference annotation to a one-transcript-per-gene canonical backbone with [`AGAT`](https://github.com/NBISweden/AGAT) + [`gffread`](https://github.com/gpertea/gffread). The canonical backbone is used by all ORF callers, P-site calibration, and DTE; the full GTF stays in scope for genome-guided alignment ([#161](https://github.com/nf-core/riboseq/issues/161)).
2. (default) Check reads distribution around annotated protein coding regions, show frame bias and estimate P-site offset for different read-length groups ([`Ribo-TISH quality`](https://github.com/zhpn1024/ribotish)).
3. (default) Predict translated ORFs and TIS de novo from alignment data ([`Ribo-TISH predict`](https://github.com/zhpn1024/ribotish)). Accepts a secondary annotation (`-a`) for background and ORF classification when extended discovery is on.
4. (default) Identify translated ORFs using P-site periodicity and read density ([`RiboCode`](https://github.com/zhengtaoxiao/RiboCode)).
5. (opt-in, `--run_ribotricer`) Derive candidate ORFs from reference data and detect translated ORFs from that list ([`Ribotricer`](https://github.com/smithlabcode/ribotricer)). Disabled by default: benchmarking found its ORF-score column is rank-unstable across biological replicates, so it is excluded from the default caller set and from cross-caller rank aggregation when enabled ([#163](https://github.com/nf-core/riboseq/issues/163)).
6. (opt-in, `--run_rpbp`) Bayesian-strict ORF calling with [Rp-Bp](https://github.com/dieterich-lab/rp-bp) (Malone et al., 2017). Score (Bayes factor) is rank-stable and participates in cross-caller rank aggregation ([#169](https://github.com/nf-core/riboseq/issues/169)).
7. (opt-in, `--run_price`) Cohort-level Bayesian ORF calling with PRICE (Erhard et al., 2018; `gedi -e Price`). Estimates a shared codon-position model across the Ribo-seq cohort by EM ([#170](https://github.com/nf-core/riboseq/issues/170)).
8. (default) Derive P-sites and QC from transcriptome alignments ([`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz)) and identify P-sites + emit bedGraph tracks ([`plastid`](https://plastid.readthedocs.io/en/latest/index.html)).
9. (opt-in, `--extended_orf_analysis true`) Discover ORFs on novel intergenic transcripts. The pipeline runs [StringTie](https://ccb.jhu.edu/software/stringtie/) reference-guided assembly (preferring RNA-seq BAMs, falling back to Ribo-seq) or accepts a user-supplied GTF (`--novel_gtf`), filters by [gffcompare](https://github.com/gpertea/gffcompare) class codes (default `u` = intergenic), and concatenates with the canonical backbone to form a hybrid annotation. Ribo-TISH predict, Ribotricer, RiboCode (via a second STAR pass against a hybrid transcriptome index), Rp-Bp and PRICE consume the hybrid annotation ([#164](https://github.com/nf-core/riboseq/issues/164), [#165](https://github.com/nf-core/riboseq/issues/165), [#171](https://github.com/nf-core/riboseq/issues/171)).
10. (opt-in, when `--extended_orf_analysis true`) Merge per-sample, per-caller ORF calls into a cohort-level ORF catalogue with class-aware collapse (transcript-ID grouping for annotated multi-exon CDS, reciprocal-overlap clustering for single-exon novel intergenic and smORFs ≤ 100 aa). Emits a unified BED12 + sidecar TSV with `called_by_<caller>` consensus columns and an AA FASTA ([#167](https://github.com/nf-core/riboseq/issues/167)).
11. (opt-in, when `--extended_orf_analysis true` + plastid) Per-ORF in-frame P-site quantification, producing an ORF × sample count matrix alongside the gene-level one ([#166](https://github.com/nf-core/riboseq/issues/166)).
12. (optional) Translational efficiency analysis at gene level with [anota2seq](https://bioconductor.org/packages/release/bioc/html/anota2seq.html) (default) or [deltaTE](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpmb.108) (DESeq2 interaction model). **Requires matched RNA-seq and Ribo-seq data**.
13. (optional, when the ORF catalogue is built) ORF-level differential translational efficiency: the deltaTE DESeq2 module is reused at ORF resolution after a join step expands the gene-level RNA-seq matrix onto ORF rows, plus a re-aggregation of canonical-CDS ORF counts into a cleaner gene-level numerator before anota2seq / deltaTE ([#168](https://github.com/nf-core/riboseq/issues/168)).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness,type
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward,riboseq
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Each row should have a 'type' value of `riboseq`, `tiseq` or `rnaseq`. Future iterations of the workflow will conduct paired analysis of matched riboseq and rnaseq samples to accomplish analysis types such as 'translational efficiency, but in the current version you should set this to `riboseq` or `tiseq` for reglar Ribo-seq or TI-seq data respectively.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/riboseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

### Including a translational efficiency analysis

![anota2seq - fold change plot](docs/images/fc.png)

In the translational efficiency analysis provided by [anota2seq](https://bioconductor.org/packages/release/bioc/html/anota2seq.html), we use matched pairs of Ribo-seq and RNA-seq data to study the relationship between transcription and translation as they differ between two treatment groups. For example the test data for this workflow has a contrasts file like:

```csv
id,variable,reference,target,batch,pair
treated_vs_control,treatment,control,treated,,pair
```

This describes how to compare groups of samples between treament groups, and between RNA-seq and Ribo-seq. In order the columns are:

- `id`: a unique identifier to use for the contrast
- 'variable`: which vaiable (column) of the sample sheet should be used to separate the treatment groups?
- `reference`: which value of the variable column should be used to select samples to be used as the reference/ base group?
- `target`: which value of the variable column should be used to select samples to be used as the target/treated group?
- `batch`: (optional) specify a variable in the sample sheet that defines sample batches
- `pair`: (optional) specify a variable in the sample sheet that defines sample pairing between RNA-seq and Ribo-seq samples. If not specified, it is assumed that the two types of sample are ordered the same.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/riboseq/usage) and the [parameter documentation](https://nf-co.re/riboseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/riboseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/riboseq/output).

## Credits

nf-core/riboseq was originally written by [Jonathan Manning](https://github.com/pinin4fjords) (Bioinformatics Engineer at Seqera) with support from [Altos Labs](https://www.altoslabs.com/) and in discussion with [Felix Krueger](https://github.com/FelixKrueger) and [Christel Krueger](https://github.com/ChristelKrueger). We thank the following people for their input:

- Anne Bresciani (ZS)
- [Felipe Almeida](https://github.com/fmalmeida) (ZS)
- [Mikhail Osipovitch](https://github.com/mosi223) (ZS)
- [Edward Wallace](https://github.com/ewallace) (University of Edinburgh)
- [Jack Tierney](https://github.com/JackCurragh) (University College Cork)
- [Maxime U Garcia](https://github.com/maxulysse) (Seqera)
- [Ira A Iosub](https://github.com/iraiosub) (The Francis Crick Institute)
- [Sebastian Uhrig](https://github.com/suhrig) (Freelance bioinformatician)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#riboseq` channel](https://nfcore.slack.com/channels/riboseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/riboseq for your analysis, please cite it using the following doi: [10.5281/zenodo.10966364](https://doi.org/10.5281/zenodo.10966364)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
