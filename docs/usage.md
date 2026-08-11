# nf-core/riboseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/riboseq/usage](https://nf-co.re/riboseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes. If you set the strandedness value to `auto` the pipeline will sub-sample the input FastQ files to 1 million reads, use Salmon Quant to infer the strandedness automatically and then propagate this information to the remainder of the pipeline. If the strandedness has been inferred or provided incorrectly a warning will be present at the top of the MultiQC report so please be sure to check when looking at the QC for your samples.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness,type
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto,riboseq
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto,riboseq
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto,riboseq
```

### Linting

By default, the pipeline will run [fq lint](https://github.com/stjude-rust-labs/fq) on all input FASTQ files, both at the start of preprocessing and after each preprocessing step that manipulates FASTQ files. If errors are found, and error will be reported and the workflow will stop.

The `extra_fqlint_args` parameter can be manipulated to disable [any validator](https://github.com/stjude-rust-labs/fq?tab=readme-ov-file#validators) from `fq` you wish. For example, we have found that checks on the names of paired reads are prone to failure, so that check is disabled by default (setting `extra_fqlint_args` to `--disable-validator P001`).

### Strandedness Prediction

If you set the strandedness value to `auto`, the pipeline will sub-sample the input FastQ files to 1 million reads, use Salmon Quant to automatically infer the strandedness, and then propagate this information through the rest of the pipeline. This behavior is controlled by the `--stranded_threshold` and `--unstranded_threshold` parameters, which are set to 0.8 and 0.1 by default, respectively. This means:

- **Forward stranded:** At least 80% of the fragments are in the 'forward' orientation.
- **Unstranded:** The forward and reverse fractions differ by less than 10%.
- **Undetermined:** Samples that do not meet either criterion, possibly indicating issues such as genomic DNA contamination.

**Note:** These thresholds apply to both the strandedness inferred from Salmon outputs for input to the pipeline and how strandedness is inferred from RSeQC results using pipeline outputs.

#### Usage Examples

1. **Forward Stranded Sample:**
   - Forward fraction: 0.85
   - Reverse fraction: 0.15
   - **Classification:** Forward stranded

2. **Reverse Stranded Sample:**
   - Forward fraction: 0.1
   - Reverse fraction: 0.9
   - **Classification:** Reverse stranded

3. **Unstranded Sample:**
   - Forward fraction: 0.45
   - Reverse fraction: 0.55
   - **Classification:** Unstranded

4. **Undetermined Sample:**
   - Forward fraction: 0.6
   - Reverse fraction: 0.4
   - **Classification:** Undetermined

You can control the stringency of this behavior with `--stranded_threshold` and `--unstranded_threshold`.

#### Errors and Reporting

The results of strandedness inference are displayed in the MultiQC report under 'Strandedness Checks'. This shows any provided strandedness and the results inferred by both Salmon (when strandedness is set to 'auto') and RSeQC. Mismatches between input strandedness (explicitly provided by the user or inferred by Salmon) and output strandedness from RSeQC are marked as fails. For example, if a user specifies 'forward' as strandedness for a library that is actually reverse stranded, this is marked as a fail.

![MultiQC - Strand check table](images/mqc_strand_check.png)

Be sure to check the strandedness report when reviewing the QC for your samples.

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 5 columns to match those defined in the table below.

A final samplesheet file consisting of both single-end Ribo-seq samples and paired-end RNA-seq data may look something like the one below.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness,type,sample_description,pair,treatment
SRX11780879,SRX11780879_SRR15480782_chr20_1.fastq.gz,SRX11780879_SRR15480782_chr20_2.fastq.gz,auto,rnaseq,PM2_5_0_1,1,control
SRX11780880,SRX11780880_SRR15480783_chr20_1.fastq.gz,SRX11780880_SRR15480783_chr20_2.fastq.gz,auto,rnaseq,PM2_5_0_2,2,control
SRX11780881,SRX11780881_SRR15480784_chr20_1.fastq.gz,SRX11780881_SRR15480784_chr20_2.fastq.gz,auto,rnaseq,PM2_5_0_3,3,control
SRX11780882,SRX11780882_SRR15480785_chr20_1.fastq.gz,SRX11780882_SRR15480785_chr20_2.fastq.gz,auto,rnaseq,PM2_5_400_1,4,treated
SRX11780883,SRX11780883_SRR15480786_chr20_1.fastq.gz,SRX11780883_SRR15480786_chr20_2.fastq.gz,auto,rnaseq,PM2_5_400_2,5,treated
SRX11780884,SRX11780884_SRR15480787_chr20_1.fastq.gz,SRX11780884_SRR15480787_chr20_2.fastq.gz,auto,rnaseq,PM2_5_400_3,6,treated
SRX11780885,SRX11780885_SRR15480788_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_C01,1,control
SRX11780886,SRX11780886_SRR15480789_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_C02,2,control
SRX11780887,SRX11780887_SRR15480790_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_C03,3,control
SRX11780888,SRX11780888_SRR15480791_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_P4001,4,treated
SRX11780889,SRX11780889_SRR15480792_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_P4002,5,treated
SRX11780890,SRX11780890_SRR15480793_chr20_1.fastq.gz,,auto,riboseq,Ribo-seq_P4003,6,treated
```

| Column         | Description                                                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandedness` | Sample strand-specificity. Must be one of `unstranded`, `forward`, `reverse` or `auto`.                                                                                                |
| `type`         | Type of sample. Must be one of `riboseq`, `rnaseq` or `tiseq`                                                                                                                          |
| `with_umi`     | (Optional) Whether the sample contains UMIs. Must be `true` or `false`. Overrides the global `--with_umi` setting for this sample.                                                     |
| `trim_length`  | (Optional) Target read length for read length equalisation. See [Read length equalisation](#read-length-equalisation).                                                                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

When specifying `contrasts` to perform a translational efficiency analysis (see below), the `type` column is necessary to distinguish Ribo-seq and RNA-seq samples. There must further be a column somewhere in the table that separates the treatment groups to be compared (`treatment` in the above example). Optionally sample pairing can be specified via an additional column (`pair` in the example), by default the same ordering will be assumed between RNA-seq and Riboseq samples of the respective groups.

## Adapter trimming options

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) performs quality and adapter trimming on FastQ files, and will automatically detect and trim the appropriate adapter sequence. The 2.x series is a self-contained Rust program with both adapter trimming and FastQC reporting built in, so it no longer wraps external Cutadapt or FastQC installations. It is the default trimming tool used by this pipeline, however you can use fastp instead by specifying the `--trimmer fastp` parameter. [fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast, all-in-one preprocessing for FastQ files. It has been developed in C++ with multithreading support to achieve higher performance. You can specify additional options for Trim Galore! and fastp via the `--extra_trimgalore_args` and `--extra_fastp_args` parameters, respectively.

With `--trimmer fastp`, `--fastp_merge` stitches overlapping paired-end mates into a single merged read. It has no effect on single-end data or with `--trimmer trimgalore`.

> **NB:** The pipeline reserves threads for Trim Galore!'s non-worker roles, passing `--cores` as `task.cpus - 3` for single-end and `task.cpus - 4` for paired-end data, clamped to between 1 and 8. Multi-core trimming therefore requires more than 4 CPUs for single-end data and more than 5 for paired-end data, and the cap of 8 worker cores is reached at 11 and 12 CPUs respectively. Trim Galore! 2.x uses an N+4 thread model (N workers plus two decompressors, a batcher and a writer) and scales near-linearly up to `--cores 8`, beyond which gzip output I/O usually becomes the limiting factor rather than worker capacity. See the [Trim Galore! changelog](https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md) and the [discussion whilst adding this logic to the nf-core/atacseq pipeline](https://github.com/nf-core/atacseq/pull/65).

## rRNA removal options

Ribosomal RNA (rRNA) removal is enabled by default (`--remove_ribo_rna`). The pipeline supports three different tools for rRNA removal, selectable via the `--ribo_removal_tool` parameter.

> [!TIP]
> For tools that use a reference database (SortMeRNA and Bowtie2), although rRNA is the primary target, the reference database can include additional abundant contaminant sequences commonly removed in ribosome profiling workflows, such as tRNAs or other non-coding RNAs. Simply add the paths to your custom FASTA files in the manifest file.

### SortMeRNA (default)

[SortMeRNA](https://github.com/biocore/sortmerna) uses k-mer matching against rRNA databases to identify and filter rRNA reads. This is the default option and requires an rRNA database manifest file.

```bash
nextflow run nf-core/riboseq --ribo_removal_tool sortmerna ...
```

By default, [rRNA databases](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) defined in the SortMeRNA GitHub repo are used. You can see an example in the pipeline GitHub repository in `assets/rrna-db-defaults.txt` which is used by default via the `--ribo_database_manifest` parameter.

> [!NOTE]
> The default databases are based on SILVA 119, which requires [licensing for commercial use](https://www.arb-silva.de/silva-license-information). SILVA 138+ uses CC-BY 4.0 licensing that freely permits commercial use with attribution. If you have licensing concerns, consider using Bowtie2 with custom rRNA reference sequences via `--ribo_removal_tool bowtie2`.

### Bowtie2

[Bowtie2](https://github.com/BenLangmead/bowtie2) performs alignment-based filtering against rRNA reference sequences. Reads that align to the rRNA references are filtered out, and unaligned reads are kept for downstream analysis. This option also requires an rRNA database manifest file specified via `--ribo_database_manifest`.

```bash
nextflow run nf-core/riboseq --ribo_removal_tool bowtie2 ...
```

### RiboDetector

> [!WARNING]
> RiboDetector has known issues with ONNX multiprocessing that can cause hangs in containerized environments (Docker, Singularity). This makes it unreliable for production use in Nextflow pipelines. We recommend using SortMeRNA or Bowtie2 for rRNA removal until these issues are resolved upstream. See [hzi-bifo/RiboDetector#61](https://github.com/hzi-bifo/RiboDetector/pull/61) for details.

[RiboDetector](https://github.com/hzi-bifo/RiboDetector) uses machine learning to identify rRNA reads without requiring a reference database. This makes it particularly useful when working with organisms that lack well-characterized rRNA sequences, or when you want to avoid database licensing requirements.

```bash
nextflow run nf-core/riboseq --ribo_removal_tool ribodetector ...
```

RiboDetector automatically determines read length from your data and uses its pre-trained neural network model to classify reads.

#### Memory management

Without memory controls, RiboDetector can consume very large amounts of RAM with large datasets. The pipeline uses the `--ribodetector_chunk_size` parameter (default: 100) to control memory usage by limiting how many reads are loaded at once. Each chunk unit represents 1024 reads, so the default of 100 loads ~102,400 reads at a time, typically resulting in ~8GB memory usage.

If you encounter memory issues, you can lower this value:

```bash
nextflow run nf-core/riboseq --ribo_removal_tool ribodetector --ribodetector_chunk_size 50 ...
```

Conversely, if you have ample memory and want faster processing, you can increase it. Set to `0` or `null` to disable chunking entirely (not recommended for large datasets).

## Read length equalisation

When comparing RNA-seq and Ribo-seq data for translational efficiency analysis, the read lengths differ substantially: RNA-seq reads are typically 75-150bp while Ribo-seq ribosome-protected fragments are 26-34bp. The pipeline provides an optional read length equalisation feature that trims RNA-seq reads to match Ribo-seq lengths before quantification. This can be enabled with the `--equalise_read_lengths` parameter.

### How it works

1. **Target length determination**: For each RNA-seq sample, the target trim length is determined by (in priority order):
   - Per-sample `trim_length` column in the samplesheet
   - Global `--equalise_read_lengths_target` parameter
   - Average length derived from the paired Ribo-seq sample (via `seqkit stats`)
2. **Ribo-seq read length measurement**: `seqkit stats` is only run on Ribo-seq samples that have paired RNA-seq samples needing length derivation (i.e., when neither a per-sample `trim_length` nor global target is specified)
3. **Trimming**: RNA-seq reads are hard-trimmed from the 5' end using TrimGalore's `--hardtrim5` option
4. **Paired-end handling**: For paired-end RNA-seq, only R1 is retained (preserving 5' position information)

If none of these sources provide a trim length for an RNA-seq sample, the pipeline will exit with an error.

### When to use

By default, the pipeline does **not** equalise read lengths. The pipeline uses STAR for alignment followed by Salmon in alignment-based mode, which applies effective length normalisation at quantification.

However, for translational efficiency (TE) analysis, read length differences can introduce bias. TE is calculated as the ratio of Ribo-seq to RNA-seq abundance - a statistical interaction, not just a comparison. When 30nt Ribo-seq reads and 150nt RNA-seq reads map to different "effective transcriptomes" (due to mappability differences), the ratio itself becomes skewed. Regions may appear translationally silent simply because short Ribo-seq reads couldn't map uniquely there.

Consider enabling `--equalise_read_lengths` if:

- You are performing TE analysis (deltaTE, anota2seq, Xtail, Riborex) and want matched mappability between modalities
- Your analysis protocol requires matched read lengths for methodological consistency
- You are replicating methods from publications that use this approach

> **Alternative approach**: Instead of trimming, you can use `--te_quantification_method pseudo` which runs Salmon pseudo-alignment for both modalities using the same k-mer index. This may help address length-related quantification biases without discarding sequence information. See [Quantification method](#quantification-method) for details.

### Trade-offs

When using read length equalisation, be aware that:

- For paired-end RNA-seq, only R1 is retained after trimming
- Shorter reads may have higher multi-mapping rates
- Some information from the original RNA-seq reads is discarded

### Example usage

```bash
nextflow run nf-core/riboseq \
    --equalise_read_lengths \
    --input samplesheet.csv \
    ...
```

Or with a global target length (useful when you don't have paired samples):

```bash
nextflow run nf-core/riboseq \
    --equalise_read_lengths \
    --equalise_read_lengths_target 28 \
    --input samplesheet.csv \
    ...
```

## Alignment options

The pipeline currently uses [STAR](https://github.com/alexdobin/STAR) to map the raw FastQ reads to the reference genome and project the alignments onto the transcriptome. STAR is fast but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

### Unique Molecular Identifiers (UMIs)

The pipeline supports UMIs to increase the accuracy of the quantification. UMIs are short sequences used to uniquely tag each molecule in a sample library and facilitate the accurate identification of read duplicates. They must be added during library preparation and prior to sequencing, therefore require appropriate arrangements with your sequencing provider.

To take UMIs into consideration for every sample in a workflow run, specify the `--with_umi` parameter. For a mixture of UMI and non-UMI libraries, add a `with_umi` column to the samplesheet and set each row to `true` or `false`. A samplesheet value overrides the global parameter for that sample. Multiple sequencing runs belonging to one sample must use the same value.

The pipeline supports UMIs embedded within a read's sequence and UMIs whose sequence is given inside the read's name. Please consult your kit's manual and/or contact your sequencing provider regarding the exact specification. UMI extraction, `--umi_discard_read`, and deduplication apply only to UMI-enabled samples. Their settings remain global across those samples, so a single run cannot combine libraries that require different barcode patterns or grouping methods.

For example, this samplesheet extracts and deduplicates the Ribo-seq library while processing the matched RNA-seq library without UMI handling:

```csv
sample,fastq_1,fastq_2,strandedness,type,with_umi
sample_ribo,ribo_R1.fastq.gz,,forward,riboseq,true
sample_rna,rna_R1.fastq.gz,rna_R2.fastq.gz,reverse,rnaseq,false
```

The samplesheet column and global parameter interact as follows:

| Configuration                                       | Result                                                                  |
| --------------------------------------------------- | ----------------------------------------------------------------------- |
| No `with_umi` column                                | Every sample uses the global `--with_umi` value.                        |
| `with_umi` is `true` or `false` for a sample        | The samplesheet value overrides the global value for that sample.       |
| `--with_umi` with selected samples set to `false`   | UMI processing applies to all samples except those explicitly disabled. |
| No `--with_umi` with selected samples set to `true` | UMI processing applies only to the samples explicitly enabled.          |

If UMIs are already embedded in the read names, mark the applicable samples with `with_umi=true` and use `--skip_umi_extract`. The pipeline will skip sequence extraction but will still deduplicate the marked samples by UMI.

The `--umitools_grouping_method` parameter affects [how similar, but non-identical UMIs](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html#method) are treated. `directional`, the default setting, is most accurate, but computationally very demanding. Consider `percentile` or `unique` if processing many samples.

#### Examples:

| UMI type     | Source                                                                                                                                                                                                                                              | Pipeline parameters                                                                                                                            |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| In read name | [Illumina BCL convert >3.7.5](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf)                                     | `--with_umi --skip_umi_extract --umitools_umi_separator ":"`                                                                                   |
| In sequence  | [Lexogen QuantSeq® 3’ mRNA-Seq V2 FWD](https://www.lexogen.com/quantseq-3mrna-sequencing) + [UMI Second Strand Synthesis Module](https://faqs.lexogen.com/faq/how-can-i-add-umis-to-my-quantseq-libraries)                                          | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern "^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"`                                   |
| In sequence  | [Lexogen CORALL® Total RNA-Seq V1](https://www.lexogen.com/corall-total-rna-seq/)<br> > _mind [Appendix H](https://www.lexogen.com/wp-content/uploads/2020/04/095UG190V0130_CORALL-Total-RNA-Seq_2020-03-31.pdf) regarding optional trimming_       | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern "^(?P<umi_1>.{12}).*"`<br>Optional: `--clip_r2 9 --three_prime_clip_r2 12` |
| In sequence  | [Takara Bio SMARTer® Stranded Total RNA-Seq Kit v3](https://www.takarabio.com/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA/SMARTer%20Stranded%20Total%20RNA-Seq%20Kit%20v3%20-%20Pico%20Input%20Mammalian%20User%20Manual-a_114949.pdf) | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern2 "^(?P<umi_1>.{8})(?P<discard_1>.{6}).*"`                                  |

> _No warranty for the accuracy or completeness of the parameters is implied_

## Reference genome options

Please refer to the [nf-core website](https://nf-co.re/usage/reference_genomes) for general usage docs and guidelines regarding reference genomes.

### Explicit reference file specification (recommended)

The minimum reference genome requirements for this pipeline are a FASTA and GTF file, all other files required to run the pipeline can be generated from these files. For example, the latest reference files for human can be derived from Ensembl like:

```
latest_release=$(curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | grep -o '"release":[0-9]*' | cut -d: -f2)
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${latest_release}.gtf.gz
```

These files can then be specified to the workflow with the `--fasta` and `--gtf` parameters.

Notes:

- Compressed reference files are supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

- If `--gff` is provided as input then this will be converted to a GTF file, or the latter will be used if both are provided.
- If `--gene_bed` is not provided then it will be generated from the GTF file.
- If `--additional_fasta` is provided then the features in this file (e.g. ERCC spike-ins) will be automatically concatenated onto both the reference FASTA file as well as the GTF annotation before building the appropriate indices.

#### Indices

By default, indices are generated dynamically by the workflow for tools such as STAR and Salmon. Since indexing is an expensive process in time and resources you should ensure that it is only done once, by retaining the indices generated from each batch of reference files:

- the `--save_reference` parameter will save your indices in your results directory

Once you have the indices from a workflow run you should save them somewhere central and reuse them in subsequent runs using custom config files or command line parameters such as `--star_index '/path/to/STAR/index/'`.

#### Gencode

If you are using [GENCODE](https://www.gencodegenes.org/) reference genome files please specify the `--gencode` parameter because the format of these files is slightly different to ENSEMBL genome files:

- The `--gtf_group_features_type` parameter will automatically be set to `gene_type` as opposed to `gene_biotype`, respectively.
- If you are running Salmon, the `--gencode` flag will also be passed to the index building step to overcome parsing issues resulting from the transcript IDs in GENCODE fasta files being separated by vertical pipes (`|`) instead of spaces (see [this issue](https://github.com/COMBINE-lab/salmon/issues/15)).

### iGenomes (not recommended)

If the `--genome` parameter is provided (e.g. `--genome GRCh37`) then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.

However this is no longer recommended because:

- Gene annotations in iGenomes are extremely out of date. This can be particularly problematic for RNA-seq analysis, which relies on accurate gene annotation.
- Some iGenomes references (e.g., GRCh38) point to annotation files that use gene symbols as the primary identifier. This can cause issues for downstream analysis, such as the nf-core [differential abundance](https://nf-co.re/differentialabundance) workflow where a conventional gene identifier distinct from symbol is expected.

### Custom genome stanzas

`--genome` is not restricted to iGenomes keys: any config file supplied with `-c` can define its own `genomes` block, and the pipeline will take the attributes below from the stanza matching `--genome`.

| Genome attribute   | Equivalent parameter |
| ------------------ | -------------------- |
| `fasta`            | `--fasta`            |
| `gtf`              | `--gtf`              |
| `gff`              | `--gff`              |
| `transcript_fasta` | `--transcript_fasta` |
| `additional_fasta` | `--additional_fasta` |
| `star`             | `--star_index`       |
| `salmon`           | `--salmon_index`     |
| `kallisto`         | `--kallisto_index`   |
| `bbsplit`          | `--bbsplit_index`    |
| `sortmerna`        | `--sortmerna_index`  |

```groovy title="my_genomes.config"
params {
    genomes {
        'GRCh38_local' {
            fasta     = '/path/to/genome.fa'
            gtf       = '/path/to/genes.gtf'
            star      = '/path/to/star/'
            salmon    = '/path/to/salmon/'
            sortmerna = '/path/to/sortmerna/'
        }
    }
}
```

```bash
nextflow run nf-core/riboseq -profile docker -c my_genomes.config --genome GRCh38_local --input samplesheet.csv --outdir results
```

Setting the equivalent parameter explicitly (on the command line, in a `-params-file`, or in a config `params` block) overrides the genome attribute.

### GTF filtering

By default, the input GTF file will be filtered to ensure that sequence names correspond to those in the genome fasta file, and to remove rows with empty transcript identifiers. Filtering can be bypassed completely where you are confident it is not necessary, using the `--skip_gtf_filter` parameter. If you just want to skip the 'transcript_id' checking component of the GTF filtering script used in the pipeline this can be disabled specifically using the `--skip_gtf_transcript_filter` parameter.

### Canonical annotation backbone

Ribo-seq footprint reads are ~28-32 nt and cannot resolve which isoform is being translated when reads fall on exons shared between isoforms (see [Wang et al. 2016, Bioinformatics 32:1880](https://academic.oup.com/bioinformatics/article/32/12/1880/1744291) for a quantitative analysis: only ~7% of human Ribo-seq reads map uniquely against a full multi-isoform annotation). The pipeline therefore separates the annotation used for genome-guided alignment (full multi-isoform, supplied via `--gtf`) from the annotation backbone used by the genome-coordinate ORF callers (Ribo-TISH, Ribotricer), plastid P-site quantification and the translational-efficiency analysis (one-transcript-per-gene, supplied via `--canonical_gtf`).

The transcriptome-coordinate tools (RiboCode, riboWaltz and Salmon) read the reference-transcriptome alignment that STAR emits with `--quantMode TranscriptomeSAM`, which is keyed to the full multi-isoform annotation; they therefore continue to use `--gtf`, since a canonical-only annotation would not match the transcript IDs in those BAMs.

`--canonical_gtf` is a recommended backbone, not a hard requirement: nothing enforces exactly one transcript per gene, so an annotation that occasionally lists more than one (e.g. MANE, which includes MANE Plus Clinical alongside MANE Select for some genes) is handled without error. Such genes simply retain a little of the shared-exon ambiguity the backbone is designed to reduce.

Recommended sources for `--canonical_gtf`:

| Organism                                              | Source                                                                                               | Extraction                                                 |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| Human (GRCh38)                                        | [MANE Select GTF](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/)                      | use directly                                               |
| Mouse, zebrafish, any Ensembl organism (release ≥104) | [Ensembl GTF](https://www.ensembl.org/info/genome/genebuild/canonical.html)                          | `grep 'tag "Ensembl_canonical"' input.gtf > canonical.gtf` |
| Any organism (fallback)                               | full `--gtf` + [AGAT](https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html) | done automatically when `--canonical_gtf` is omitted       |

If `--canonical_gtf` is not supplied the pipeline runs [`agat_sp_keep_longest_isoform.pl`](https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html) on the full GTF and uses the result as the backbone. AGAT operates structurally (it keeps the isoform with the longest CDS per gene, falling back to the longest concatenated exons for genes with no CDS) rather than from curation, so a curated source is strongly preferred where available.

MANE Select vs `Ensembl_canonical` for non-coding genes: MANE Select covers virtually all protein-coding genes plus a growing but partial set of non-coding genes ([NCBI Insights, MANE v1.4](https://ncbiinsights.ncbi.nlm.nih.gov/2024/10/28/mane-v1-4-mane-select-non-coding-genes/)). `Ensembl_canonical` has broader biotype priority including lncRNA and other ncRNA biotypes, so users working primarily on smORFs in lncRNAs may get better transcript recall from `Ensembl_canonical`.

## Riboseq-specific options

The pipeline will by default run the [Ribo-TISH](https://github.com/zhpn1024/ribotish) [quality](https://github.com/zhpn1024/ribotish?tab=readme-ov-file#quality) and [predict](https://github.com/zhpn1024/ribotish?tab=readme-ov-file#predict) commands for QC and ORF prediction, respectively. Additional arguments can be supplied to either command via the `--extra_ribotish_quality_args` and `--extra_ribotish_predict_args` parameters.

Ribo-TISH `quality` is fed the canonical annotation (`--canonical_gtf`, or the AGAT-derived longest-isoform fallback) rather than the full multi-isoform GTF. `quality` estimates P-site offsets and read-length QC against CDS-bearing canonical transcripts; mixing in CDS-absent or near-duplicate isoforms degrades that calibration without adding diagnostic signal. The `predict` step receives the same canonical input by default; the wrapper also accepts an optional secondary annotation on `-a` for novel-transcript discovery modes.

### ORF calling and cross-caller agreement

By default the pipeline calls ORFs with two tools, Ribo-TISH `predict` and RiboCode, and reports an ORF as agreed only when both callers support it. This intersection is precision-weighted: the agreed set is conservative and may omit ORFs that only one caller detects.

Ribotricer is available as a third caller but is off by default. Enable it with `--run_ribotricer true` for broader recall, after which an ORF is agreed on a majority vote (2 of 3). It is opt-in because its ORF-score column is unstable across biological replicates even though its binary call set is reproducible. When enabled, its binary calls count toward agreement but its score is excluded from cross-caller rank aggregation, and the pipeline warns at runtime.

### Rp-Bp (opt-in, overnight)

[Rp-Bp](https://github.com/dieterich-lab/rp-bp) (Malone et al., 2017) is a Bayesian-strict ORF caller that complements RiboCode's permissive canonical-CDS calls. It is the recommended second caller when statistical rigour matters more than turnaround time. Activate with `--run_rpbp true`.

> :warning: **Runtime cost.** Rp-Bp's Bayesian MCMC fit dominates wall-clock and takes roughly **20-24 hours per replicate** at genome-wide scale. The pipeline emits a runtime warning when `--run_rpbp` is set. Plan compute time, queue limits and instance lifetimes accordingly.

Rp-Bp's score column (Bayes factor) is stable and is retained in the cross-caller rank-aggregation set alongside RiboCode and Ribo-TISH; Ribotricer's score column is excluded due to known instability but Rp-Bp's is not.

Rp-Bp runs through the upstream `nf-core/rpbp/*` modules driven by the `FASTA_GTF_BAM_RPBP` nf-core subworkflow, which orchestrates `prepare-rpbp-genome`, `extract-metagene-profiles`, `estimate-metagene-profile-bayes-factors`, `select-periodic-offsets`, `get-periodic-lengths-offsets`, `extract-orf-profiles`, `estimate-orf-bayes-factors` and `select-final-prediction-set` from your `--fasta` / `--gtf` inputs without you having to author a YAML config. Tool CLI overrides are exposed via `--extra_rpbp_preparegenome_args` and `--extra_rpbp_predictorfs_args`.

Per-sample final-prediction outputs - filtered BED of predicted ORFs (with Bayes factor in column 5), plus matched nucleotide and protein FASTAs - are published under `<outdir>/orf_predictions/rpbp/`.

**Annotation.** Rp-Bp is given the full multi-isoform `--gtf` annotation, not the one-transcript-per-gene canonical backbone that the pipeline uses elsewhere to disambiguate P-site quantification. Rp-Bp enumerates candidate ORFs across every transcript isoform (deduplicating identical ORFs by genomic coordinate) and then resolves redundant and overlapping ORFs itself - the longest ORF per stop codon, then the highest Bayes factor among overlaps. Collapsing the annotation to one isoform per gene would silently remove ORFs that exist only on non-canonical isoforms (alternative-5'UTR uORFs, isoform-specific N-terminal extensions or truncations, retained-intron and alternative-exon ORFs) and bias the reported ORF types toward canonical CDS, with no compensating benefit; PRICE is handled the same way and for the same reason ([Malone et al., 2017](https://academic.oup.com/nar/article/45/6/2960/2953491)). Under `--extended_orf_analysis true` Rp-Bp receives the full reference with the novel intergenic genes appended, so novel transcripts come into discovery scope without giving up the isoforms.

> :information_source: **STAR alignment params vs upstream rpbp.** rpbp's own pipeline runs STAR with Ribo-seq-tuned settings (`outFilterMismatchNmax 1`, `outFilterMismatchNoverLmax 0.04`, `outFilterType BySJout`, `sjdbOverhang 33`, `winAnchorMultimapNmax 100`, `seedSearchStartLmaxOverLread 0.5`). We use the pipeline's standard STAR alignment (shared with the RNA-seq side of paired runs), which is more permissive. Practical impact: rpbp processes whatever alignments it gets, but periodicity / Bayes-factor distributions will differ from a standalone rpbp run on the same FASTQs. If you need bit-identical-to-standalone-rpbp output, override with `--extra_star_align_args '--outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterType BySJout --winAnchorMultimapNmax 100 --seedSearchStartLmaxOverLread 0.5'`. Note that `sjdbOverhang` is baked into the STAR index and cannot be changed post-hoc - it would require regenerating the index with `--sjdbOverhang 33`, and that change would only be appropriate for a Ribo-seq-only run (RNA-seq reads are too long for that setting). Tracked for future work: [#173](https://github.com/nf-core/riboseq/issues/173).

### PRICE (opt-in)

[PRICE](https://github.com/erhard-lab/gedi/wiki/Price) (Erhard et al., 2018) is a Bayesian ORF caller distributed as part of the [Gedi](https://github.com/erhard-lab/gedi) Java framework. Unlike the per-sample callers, PRICE estimates a shared codon-position model across the riboseq cohort by EM and is invoked one-shot rather than per-sample. Activate with `--run_price true`.

> :warning: **Runtime cost.** PRICE's EM fit is comparable in wall-clock to other heavy ORF callers at genome-wide scale. The pipeline emits a runtime warning when `--run_price` is set. Plan compute accordingly.

The pipeline builds a binary `.oml` genome index via `gedi -e IndexGenome` once per run, then calls PRICE once across the cohort with the index plus the riboseq BAMs. PRICE's primary output is `${prefix}.orfs.tsv`, a table of all called ORFs with start-codon score, range score, p-value (uncorrected) and per-condition / total read counts. Tool CLI arguments can be appended via `--extra_price_indexgenome_args` and `--extra_price_price_args`.

**Annotation.** Like Rp-Bp, PRICE is given the full multi-isoform `--gtf` annotation rather than the one-transcript-per-gene canonical backbone: it resolves overlapping ORFs and rescues multimappers with its own EM, so restricting it to a single isoform per gene would only narrow ORF discovery and bias ORF-type classification toward canonical CDS.

When `--extended_orf_analysis true` is set, PRICE's IndexGenome is built from the full reference with the novel intergenic genes appended, so ORFs on novel transcripts come into scope without giving up the isoforms.

PRICE's CLI banner reports `Price version 1.0.4` while the Bioconda package is `gedi 1.0.6a` (Price is one tool inside the Gedi umbrella). The pipeline captures the package version via `gedi -e Version` for `versions.yml`.

> :information_source: **Data scale.** PRICE's ORF inference step (`PriceOrfInference`) requires more candidate ORFs than the chr20-only `-profile test` data provides; on small test inputs it fails late in the pipeline with `Index out of bounds`. End-to-end PRICE validation needs realistic Ribo-seq depth and is exercised on the full-scale Platform iteration, not in the chr20 CI test set.

## P-site identification

The pipeline will by default run [riboWaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz) for P-site identification and diagnostics, unless disabled with `--skip_ribowaltz`. Additional arguments can be supplied via `--extra_ribowaltz_args` parameters. An example is: `--extra_ribowaltz_args "--length_range 27:31 --periodicity_threshold 40 --extremity 5end --start_nts 45 --stop_nts 24"`. If not provided, defaults used in the [nf-core module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ribowaltz/templates/ribowaltz.r) are used.

In addition, the pipeline runs [plastid](https://plastid.readthedocs.io/en/latest/index.html) to identify P-sites and create bedGraph tracks, unless disabled with `--skip_plastid`.

## Translational efficiency

If you have paired RNA-seq and Riboseq samples, you can use this workflow to initiate a translational efficiency analysis.

The pipeline supports two methods for translational efficiency analysis:

### anota2seq (default)

Translational efficiency analysis as conducted by [anota2seq](https://bioconductor.org/packages/release/bioc/html/anota2seq.html) involves the integrated analysis of RNA-seq and Ribo-seq data to discern changes in translational efficiency across different experimental conditions. It quantitatively assesses how variations in mRNA abundance and ribosome occupancy lead to alterations in protein synthesis, enabling the identification of genes with post-transcriptional and translational regulation.

### deltaTE

Alternatively, you can use the deltaTE method by specifying `--translational_efficiency_method deltate`. The deltaTE method, based on [Chothani et al. (2019)](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpmb.108), uses DESeq2 with an interaction model to detect differentially translated genes (DTEGs). It integrates Ribo-seq and RNA-seq data to identify genes with significant changes in translational efficiency, classifying them into biological categories based on their regulatory patterns.

Both methods analyze differences between conditions for RNA-seq and Ribo-seq samples, but use different statistical frameworks to identify translational regulation.

### ORF-level differential translation

When `--extended_orf_analysis true` is set with `--te_quantification_method plastid_psite` (the default), at least one ORF caller is enabled, and `--contrasts` is supplied, the pipeline runs differential translation at ORF resolution, additively to the gene-level DTE above and using the same method selector (`--translational_efficiency_method`). Three methods are available:

- `anota2seq` (default, gene-level + ORF-level): APV + RVM regression of Ribo on RNA per identifier.
- `deltate`: DESeq2 with a `~ condition + seq_type + condition:seq_type` interaction model.
- `dotseq` (ORF-level only): DOTSeq's ORF-level differential translation efficiency (DESeq2 + ashr) AND its DOTSeq-specific DOU contrast - a per-gene beta-binomial GLM modelling whether each ORF gains or loses a share of its parent gene's ribosome occupancy across conditions (the question "is this ORF gaining ribosomes at the expense of its siblings?", which DTE alone can't answer). Selecting `dotseq` requires `--extended_orf_analysis true` and an enabled ORF caller; the gene-level fit is skipped.

A pre-processing step joins the per-ORF P-site count matrix (`<outdir>/quantification/orf_level/orf_psite_counts.tsv`) with a gene-level RNA-seq matrix via `orf_to_gene.tsv` from the catalogue, producing one combined count table whose rows are ORFs and whose RNA columns hold the host gene's count replicated across all ORFs sharing that gene. The selected method is then fitted, treating each ORF row as a feature. Results land under `<outdir>/translational_efficiency/orf_level/<method>/`, alongside the shared combined input matrix (`<outdir>/translational_efficiency/orf_level/orf_combined_counts.tsv`).

Which RNA-seq matrix supplies the denominator depends on whether a novel-transcript source is configured.

With one (`--skip_stringtie false` or `--novel_gtf`), the RNA-seq reads are requantified against the full reference transcriptome augmented with the novel intergenic transcripts, using the same STAR-alignment then Salmon path as the primary quantification. This is what lets ORFs on novel genes participate: their host gene does not exist in the canonical reference, so a canonical-only denominator would give them no RNA row and the join would silently drop every such ORF, leaving only novel ORFs that happen to sit on known genes. Canonical genes keep their full-isoform counts (the one-transcript-per-gene backbone used for ORF calling would undercount multi-isoform genes). The augmented matrices are published under `<outdir>/quantification/salmon_hybrid/`. Novel-transcript ORFs are still subject to the low-count caveat below, and their reliability depends on the confidence of the underlying StringTie assembly.

Without one there are no novel genes, so the primary Salmon matrix (`<outdir>/quantification/salmon/`) already covers every host gene and is reused directly rather than a second STAR index being built. The ORF-level denominator is then exactly the gene-level one. It carries the Ribo-seq columns too; the join drops them.

Two caveats apply:

- **Row independence.** Multiple ORFs from the same gene share a single gene-level RNA-seq denominator row (sibling ORFs sit on the same mRNA, so there is no separately measurable per-ORF RNA level). After the join, those rows are perfectly correlated, and both anota2seq's per-identifier APV regression and deltaTE's per-row DESeq2 fit treat each ORF as an independent observation. Treat p-values for ORFs sharing a host gene as a ranking, not strict significance.
- **Low-count ORFs.** uORFs, smORFs and low-abundance novel ORFs frequently have sparse P-site counts. DESeq2 dispersion estimation and anota2seq's RVM are both unreliable on rows with too many zeros. Pass method-side options via `--extra_orf_deltate_args` (deltaTE), `--extra_orf_anota2seq_run_args` (anota2seq) or `--extra_dotseq_args` (DOTSeq); inspect the diagnostic plots before interpreting results.
- **Multiple-testing scale.** Per-ORF runs typically have 5-50x more rows than per-gene runs, so BH-adjusted p-values are more conservative by construction. This is expected and applies identically to both methods.

### Method comparison

**anota2seq** studies differences between conditions for both RNA-seq and Ribo-seq samples. It also assesses combined results from two measures as they relate to one another:

- **mRNA abundance**: Changes in total RNA levels that lead to corresponding changes in translation
- **Translation**: Differences in translation not occurring as a result of overall RNA levels
- **Buffering**: Changes in total RNA levels that do not lead to increased translation

**deltaTE** classifies genes based on statistical significance patterns:

- **mRNA_abundance**: RNA changes forwarded to translation without net translational efficiency changes
- **Translation**: Pure translational regulation - ribosome changes without mRNA changes
- **Buffering**: Translation dampens RNA changes (opposite directional effects)
- **Intensified**: Translation amplifies RNA changes (same direction, deltaTE-specific)

This table summarizes the conceptual framework:

| Category       | RNA-seq   | Ribo-seq        | Translational Efficiency |
| -------------- | --------- | --------------- | ------------------------ |
| mRNA abundance | Changed   | Changed         | Unchanged                |
| Translation    | Unchanged | Changed         | Changed                  |
| Buffering      | Changed   | Stable/Opposite | Changed                  |
| Intensified\*  | Changed   | Amplified       | Changed                  |

\*Intensified is specific to the deltaTE method.

### Method selection

By default, the pipeline uses anota2seq for translational efficiency analysis. To use the deltaTE method instead, specify:

```bash
--translational_efficiency_method deltate
```

Both methods require the same input format and contrasts specification, but produce different output files and use different statistical approaches.

### Quantification method

The pipeline offers three methods for quantifying gene expression for TE analysis, controlled by the `--te_quantification_method` parameter:

#### In-frame P-sites (default)

```bash
--te_quantification_method plastid_psite
```

This is the default method. RNA-seq reads are quantified using the alignment-based method (STAR + Salmon). Ribo-seq reads are counted only if they map to coding regions and their predicted P-sites (as determined by Plastid) coincide with the annotated reading frame. Specifically, a Ribo-seq read is counted only when its inferred P-site aligns with the reading frame of a coding sequence (CDS) defined in the provided annotation (GTF) file.

This is the scientifically appropriate quantity for Ribo-seq. Salmon's underlying model assumes uniform coverage and a length-dependent fragment distribution and relies on long or paired-end reads to resolve multi-mapping. Ribo-seq footprints (~28-32 nt) violate all three assumptions: coverage non-uniformity is the biological signal, fragment lengths are physically constrained by the ribosome, and short footprints cannot disambiguate shared exons (see [Ribomap, Bioinformatics 2016](https://academic.oup.com/bioinformatics/article/32/12/1880/1744291)).

Caveats:

- The methods used to quantify RNA-seq reads and Ribo-seq reads are not the same, giving rise to different technical biases in the counts.
- The quantification of in-frame P-sites is subject to the periodicity of the Ribo-seq experiment. Before comparing samples against each other, it should be confirmed that the periodicities of the samples are similar.
- Quantification of overlapping ORFs is still flawed and cannot be fully deconvoluted, since the counts from one ORF can leak into the other ORF unless the periodicity efficiency is perfect (which it usually is not).
- Counts are summed per gene; ORF-level quantification is not yet supported in this mode.

#### Alignment-based

```bash
--te_quantification_method alignment
```

Reads are aligned with STAR, and Salmon quantifies from the transcriptome BAM in alignment-based mode. This was the default in earlier pipeline versions and is retained for backward compatibility and for downstream tools that expect Salmon-format transcript-level counts. Choose this if you need to reproduce results from a prior run, compare against a Salmon-based cohort, or feed counts into a tool that specifically expects Salmon output. Note that switching from `alignment` to `plastid_psite` (or vice versa) changes the actual per-gene count values - downstream comparisons across pipeline versions are not apples-to-apples if the default changed between runs.

#### Pseudo-alignment

```bash
--te_quantification_method pseudo
```

Quantifies reads directly against the transcriptome for both Ribo-seq and RNA-seq samples, without going via the genome alignment. This is an **experimental alternative** that:

- Applies the same k-mer index and quantification algorithm to both modalities
- May help reduce length-related quantification biases without requiring read trimming
- Uses a k-mer size of 23 by default, suitable for short Ribo-seq reads. Adjust with `--pseudo_aligner_kmer_size` if needed.

Consider this option when:

- You want to avoid the information loss from read length equalisation trimming
- You prefer k-mer-based quantification for methodological consistency between modalities

> **Note**: The pseudo-alignment pathway runs **in addition to** the standard STAR alignment, which is still needed for position-dependent analyses (P-sites, ribosome periodicity, ORF detection). The pseudo-alignment counts are only used for TE analysis.

##### Choosing the pseudo-aligner

`--pseudo_aligner` selects the tool used for this pathway, either `salmon` (default) or `kallisto`:

```bash
--te_quantification_method pseudo --pseudo_aligner kallisto \
    --kallisto_quant_fraglen 30 --kallisto_quant_fraglen_sd 5
```

What `--pseudo_aligner` controls:

- The tool that produces the per-sample transcript quantifications feeding the TE analysis, and the `quant_type` handed to tximport
- Which index is built for the pathway: a Salmon index (supply a pre-built one with `--salmon_index`), or a kallisto index (supply a pre-built index file with `--kallisto_index`)
- The output subdirectory, `quantification/salmon_te_pseudo` or `quantification/kallisto_te_pseudo`

What it does **not** control:

- The primary STAR genome alignment, which every position-dependent analysis depends on
- The alignment-mode Salmon quantification of the STAR transcriptome BAM, which always uses Salmon
- Strandedness inference, which always subsamples and runs Salmon and therefore always needs a Salmon index

kallisto-specific caveats:

- kallisto cannot estimate a fragment length distribution from single-end reads, and Ribo-seq libraries are single-end. `--kallisto_quant_fraglen` and `--kallisto_quant_fraglen_sd` are therefore **required**, and the pipeline exits at startup if they are missing. For Ribo-seq the fragment is the footprint itself, so the values should reflect your observed footprint length distribution (commonly around 30 nt), not a typical RNA-seq insert size. Paired-end RNA-seq samples in the same run ignore them, since kallisto infers the distribution from the pairs.
- The same two values are applied to every single-end sample in the run, so a cohort mixing libraries with very different insert size distributions is better quantified with Salmon.
- MultiQC reports kallisto's run logs rather than the richer Salmon quantification metrics.
- `--pseudo_aligner_kmer_size` is passed to `kallisto index -k`, which requires an odd value no greater than 31.

### Contrasts specification

To carry out this analysis, the pipeline must be supplied with one or more 'contrasts' describing the comparison to be made.

For example the test data for this workflow has a contrasts file like:

```csv
id,variable,reference,target,batch,pair
treated_vs_control,treatment,control,treated,,pair
```

This describes how to compare groups of samples between treatment groups, and between RNA-seq and Ribo-seq. In order the columns are:

- `id`: a unique identifier to use for the contrast
- `variable`: which variable (column) of the sample sheet should be used to separate the treatment groups?
- `reference`: which value of the variable column should be used to select samples to be used as the reference/base group?
- `target`: which value of the variable column should be used to select samples to be used as the target/treated group?
- `batch`: (optional) specify a variable in the sample sheet that defines sample batches for batch effect correction in anota2seq
- `pair`: (optional) specify a variable in the sample sheet that defines sample pairing between RNA-seq and Ribo-seq samples. If not specified, it is assumed that the two types of sample are ordered the same.

> [!NOTE]
> The analysis automatically subsets the count data to only the samples involved in each contrast. Additional anota2seq options can be passed via `--extra_anota2seq_run_args` (see parameter documentation for details).

## Novel transcript discovery (StringTie / user-supplied GTF)

The pipeline can extend the canonical reference annotation with novel intergenic transcripts, either by running [StringTie](https://ccb.jhu.edu/software/stringtie/) reference-guided assembly or by accepting a user-supplied GTF via `--novel_gtf`. The two sources feed the same downstream filtering chain and produce a hybrid annotation at `<outdir>/transcript_assembly/stringtie/hybrid_reference.gtf` (canonical backbone + filtered novel transcripts, sorted by genomic position).

### Source 1: StringTie assembly

Set `--skip_stringtie false` to enable assembly. StringTie requires RNA-seq BAMs; if your samplesheet contains only Ribo-seq samples, use `--novel_gtf` with a GTF produced from a separate RNA-seq run instead.

Per-sample StringTie GTFs are merged with `stringtie --merge` into a unified annotation.

### Source 2: User-supplied GTF (`--novel_gtf`)

Pass `--novel_gtf path/to/curated.gtf` to skip StringTie entirely and feed your own annotation through the filtering chain. Useful when you already have a tissue/cell-type-specific annotation from long-read RNA-seq or a published ORF catalogue. Setting `--novel_gtf` takes precedence over `--skip_stringtie`.

### Filtering: gffcompare class codes

The novel GTF is classified against the full reference with [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) and filtered to entries whose `class_code` is in `--gffcompare_class_codes` (default `"u"`, intergenic only).

Stranded users (stranded RNA-seq and stranded Ribo-seq libraries) can extend the filter to `--gffcompare_class_codes "u,x"` to recover translated antisense transcripts (class `x` = antisense overlap with a known locus). The default remains `u`-only because antisense classification is unreliable for non-stranded or partially-stranded protocols.

### Optional rRNA / repeat blacklist

Supply `--rrna_blacklist path/to/blacklist.bed` to drop novel transcripts overlapping known rRNA or repeat regions on the same strand (`bedtools intersect -v -s`). UCSC publishes ready-made blacklists for hg38 and mm10; for other organisms you must construct one yourself (or omit the parameter to skip the step silently).

### Knobs

- `--skip_stringtie` - default `true`. Set to `false` to run StringTie assembly (requires RNA-seq BAMs in the samplesheet).
- `--novel_gtf` - user-supplied novel GTF (bypasses StringTie when set).
- `--extra_stringtie_args` - extra args passed to per-sample StringTie.
- `--extra_stringtie_merge_args` - extra args passed to `stringtie --merge` (e.g. `'-T 1 -f 0.1'` for stricter TPM and isoform-fraction cutoffs); see the [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) for the full set of merge flags.
- `--gffcompare_class_codes` - comma-separated gffcompare class codes to retain (default `u`).
- `--rrna_blacklist` - optional BED of rRNA / repeat regions to exclude.

The hybrid GTF is exposed as a workflow channel (`hybrid_gtf` emit) and is wired into the genome-BAM ORF callers (Ribo-TISH `predict`, Ribotricer) when `--extended_orf_analysis true` is set; see [Extended ORF discovery](#extended-orf-discovery) below.

## Extended ORF discovery

By default, all ORF callers run against the canonical backbone GTF so the pipeline produces well-characterised annotated-ORF calls. To discover novel ORFs in the novel intergenic transcripts produced by StringTie or supplied via `--novel_gtf`, set `--extended_orf_analysis true`. This routes the hybrid GTF (`<outdir>/transcript_assembly/stringtie/hybrid_reference.gtf`) into the ORF callers:

- **Ribo-TISH `predict`**: hybrid GTF on `-g` (discovery target); canonical backbone on `-a` (background and ORF classification).
- **Ribotricer `prepare-orfs`**: hybrid GTF directly (Ribotricer has no secondary-annotation concept; CDS-absent novel transcripts are auto-labelled `novel`).
- **RiboCode**: hybrid GTF as the annotation source plus a hybrid transcriptome BAM produced by a second STAR alignment pass (see below).

A novel-transcript source — `--skip_stringtie false` or `--novel_gtf <path>` — is required to route hybrid annotation into the callers. If `--extended_orf_analysis true` is set without one, the pipeline warns and the callers fall back to their usual annotation, but the cross-caller ORF catalogue, per-ORF P-site quantification and ORF-level DTE are still produced, built against the full reference annotation. Note `--skip_stringtie` defaults to `true`.

Two caveats when comparing runs. Catalogue `orf_id` values are assigned by coordinate order, so adding novel transcripts renumbers every ORF downstream of them — compare runs on genomic coordinates, never on `orf_id`. And with a single enabled caller the consensus view carries no cross-caller information; with `--orf_min_callers 2` it is empty, and the pipeline warns.

Without a novel source the callers do not share one annotation: Ribo-TISH and Ribotricer read the one-transcript-per-gene canonical backbone while Rp-Bp, PRICE and RiboCode read the multi-isoform `--gtf`. An ORF that exists only on a non-canonical isoform is therefore invisible to the first group by construction, so it can never reach two-caller agreement. Keep that in mind when raising `--orf_min_callers` above 1: the threshold penalises isoform-specific ORFs rather than filtering them on evidence.

### Second STAR pass for RiboCode

RiboCode (and STAR `--quantMode TranscriptomeSAM` in general) needs a transcriptome-coordinate BAM keyed to whichever transcriptome FASTA was supplied at alignment time. The primary STAR pass is built against the reference transcriptome, so novel StringTie transcripts are invisible to it. When `--extended_orf_analysis true` is set, the pipeline therefore:

1. Extracts spliced transcript sequences from the hybrid GTF with `gffread -w` to produce a hybrid transcriptome FASTA (canonical + novel).
2. Rebuilds a STAR index against the original genome FASTA, using the hybrid GTF as `--sjdbGTFfile`.
3. Re-aligns the Ribo-seq reads against that hybrid index to obtain a hybrid transcriptome BAM, which is then fed to RiboCode in place of the reference transcriptome BAM.

The second pass roughly doubles STAR alignment compute for Ribo-seq samples and consumes additional disk for the hybrid index and BAMs. It runs only when `--extended_orf_analysis true` and a novel-transcript source are configured, and is restricted to Ribo-seq samples (RNA-seq and TI-seq do not feed RiboCode). The hybrid transcriptome FASTA and hybrid STAR index are each built once per pipeline run. Outputs are published under `<outdir>/alignment/hybrid_star/` (and `<outdir>/genome/index/hybrid_star/` for the hybrid STAR index, when `--save_reference` is set).

**riboWaltz stays on the primary reference-transcriptome BAM by design.** riboWaltz is a QC/calibration tool and its CDS-dependent plots (frame distribution, start/stop metaprofiles) are driven by annotated CDS-bearing transcripts. Routing CDS-absent novel transcripts through riboWaltz would dilute its diagnostic plots without contributing to ORF discovery (riboWaltz does not call ORFs). Salmon likewise stays on the primary reference transcriptome, and plastid P-site quantification on the canonical backbone, regardless of `--extended_orf_analysis`.

**Ribo-TISH `quality` stays on the canonical backbone.** The `quality` step estimates P-site offsets and read-length QC against CDS-bearing canonical transcripts. Feeding it the hybrid GTF would mix in CDS-absent novel transcripts that the P-site model cannot interpret. The Ribo-TISH `predict` step (which actually calls ORFs) is the one that consumes the hybrid GTF on `-g` plus the canonical backbone on `-a`. Empirical confirmation that canonical-only is the correct choice for `quality` is pending full-scale validation; the default reflects the spec recommendation.

The default `--extended_orf_analysis false` keeps the pre-#165 behaviour unchanged: the primary STAR pass is the only alignment and every ORF caller runs against the canonical backbone.

### Cross-sample ORF catalogue

When `--extended_orf_analysis true` is set and at least one ORF caller is enabled, the pipeline produces a cohort-level ORF catalogue under `<outdir>/orf_predictions/catalogue/`. The catalogue normalises each per-sample, per-caller output into a unified BED12 (genomic blocks, multi-exon-aware), then merges across samples and callers by clustering strategy. `orf_class` selects the strategy but never enters a grouping key, so callers that disagree on an ORF's class still merge into one row:

- annotated CDS are grouped by `transcript_id` to preserve intron-chain identity, then clustered by 80% reciprocal overlap within the transcript so a short truncated variant stays distinct from the full-length CDS;
- transcript-anchored ORFs (`uORF`, `uoORF`, `dORF`, `doORF`, `intORF`, `other`) are collapsed by transcript, strand and outer span;
- novel intergenic ORFs are clustered by 80% reciprocal overlap on the outer genomic span;
- small ORFs (`is_smorf`, i.e. ≤ `--smorf_max_aa` aa, default 100) are additionally peptide-level deduplicated: the catalogue amino-acid FASTA is clustered with MMseqs2 (`--min-seq-id 0.9 -c 0.8`) and each multi-member cluster is folded to one representative, following the GENCODE Ribo-seq ORF catalogue convention (Mudge et al. 2022). Pass `--skip_orf_collapse` to publish the coordinate-merged catalogue without this sequence-level collapse;
- cross-caller consensus is recorded in `called_by_<caller>` binary columns plus `score_<caller>` columns for Ribo-TISH / RiboCode / Rp-Bp / PRICE (Ribotricer scores are excluded from rank aggregation).

`orf_class` records only where an ORF sits relative to the annotated CDS, never how long it is; `is_smorf` carries the length flag. Raise `--smorf_max_aa` (e.g. `150`) to widen the microprotein window — this flags more ORFs as small and brings more of them into the peptide collapse, but leaves every `orf_class` value untouched. Not every caller can report every class: Ribo-TISH's `5'UTR` covers both `uORF` and `uoORF`, PRICE has no `doORF`, and ribotricer's `internal` is a fall-through rather than a frame-tested call so it stays `other`. `orf_type_native` preserves each caller's own label for auditing.

Host gene and transcript ids are resolved against the union of the full multi-isoform annotation (`--gtf`) and the filtered novel transcripts. The callers do not all receive the same annotation (PRICE takes the full multi-isoform reference, the genome-BAM callers the canonical backbone), so only that union covers every transcript an ORF can be reported on; an ORF called on an isoform absent from the annotation falls back to whatever gene label its caller emitted, which for PRICE is every gene its genomic span overlaps, concatenated. With a novel-transcript source the union is also the reference the ORF-level DTE RNA denominator is quantified against, so catalogue gene ids share a namespace with that matrix. Without one the denominator is the primary Salmon matrix, keyed on `--gtf` alone: a `--canonical_gtf` whose gene ids are absent from `--gtf` puts ORFs called against it outside that namespace, and the join drops them (the count is reported on `DTE_COUNTS_PREP`'s stderr).

The catalogue runs once per pipeline invocation (cohort-level, not per sample) and gates on `--extended_orf_analysis true` plus a non-empty enabled-caller set. The default-off path keeps the pre-#167 behaviour unchanged.

Alongside the full catalogue, a consensus view is published under `<outdir>/orf_predictions/catalogue/consensus/` containing only ORFs supported by at least `--orf_min_callers` distinct callers and recurring in at least `--orf_min_samples` samples (both default 1, so the consensus view equals the full catalogue out of the box). Raise either threshold (e.g. `--orf_min_callers 2`) for a higher-confidence catalogue that tames downstream ORF-level multiple testing; the full unfiltered catalogue is always published regardless. The threshold is applied after the small-ORF collapse, so the consensus is the high-confidence subset of the de-redundified catalogue and a micropeptide folded across several loci is judged on its combined cross-caller / cross-sample evidence (when `--skip_orf_collapse` is set it comes from the merged catalogue instead).

See [ORF catalogue (cross-sample)](output.md#orf-catalogue-cross-sample) in the output docs for the full list of published files.

### Per-ORF P-site quantification

When the cohort catalogue is built and plastid is enabled (`--skip_plastid false`, the default), the pipeline also produces a per-ORF in-frame P-site count matrix at `<outdir>/quantification/orf_level/orf_psite_counts.tsv`. This is an ORF x sample matrix of raw integer counts, complementing the gene-level matrix at `<outdir>/quantification/inframe_psite/gene_counts.tsv`. Frames for catalogue ORFs are defined by each ORF's own start codon (ATG = frame 0), not by the GTF `phase` field, so novel transcripts and non-canonical starts are handled correctly. The matrix is the input for the per-ORF translational-efficiency analysis tracked in #168; gene-level DTE is unchanged by this addition.

When `--skip_plastid true` is set together with `--extended_orf_analysis true`, the catalogue is still built but the per-ORF count matrix is skipped (the plastid wiggle tracks are not available), and a runtime warning is emitted.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run \
    nf-core/riboseq \
    --input <SAMPLESHEET> \
    --outdir <OUTDIR> \
    --gtf <GTF> \
    --fasta <GENOME FASTA> \
    -profile docker
```

> **NB:** Loading iGenomes configuration remains the default for reasons of consistency with other workflows, but should be disabled when not using iGenomes, applying the recommended usage above.

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/riboseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/riboseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/riboseq releases page](https://github.com/nf-core/riboseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
