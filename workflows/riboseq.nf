/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.transcript_fasta = WorkflowMain.getGenomeAttribute(params, 'transcript_fasta')
params.additional_fasta = WorkflowMain.getGenomeAttribute(params, 'additional_fasta')
params.gtf              = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff              = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed         = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.bbsplit_index    = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.star_index       = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index     = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index       = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index     = WorkflowMain.getGenomeAttribute(params, 'salmon')
params.kallisto_index   = WorkflowMain.getGenomeAttribute(params, 'kallisto')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowRiboseq.initialise(params, log)

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Determine whether to filter the GTF or not
def filterGtf =
    ((
        // Condition 1: Alignment is required and aligner is set
        !params.skip_alignment && params.aligner
    ) ||
    (
        // Condition 2: Pseudoalignment is required and pseudoaligner is set
        !params.skip_pseudo_alignment && params.pseudo_aligner
    ) ||
    (
        // Condition 3: Transcript FASTA file is not provided
        !params.transcript_fasta
    )) &&
    (
        // Condition 4: --skip_gtf_filter is not provided
        !params.skip_gtf_filter
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { PREPARE_GENOME                    } from '../subworkflows/local/prepare_genome'
include { PREPROCESS_RNASEQ                 } from '../subworkflows/nf-core/preprocess_rnaseq'
include { FASTQ_ALIGN_STAR                  } from '../subworkflows/nf-core/fastq_align_star'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                                              } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                          } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_SORT                                        } from '../modules/nf-core/samtools/sort'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/nf-core/umitools/prepareforrsem'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RIBOSEQ {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.gencode,
        is_aws_igenome,
        biotype,
        prepareToolIndices,
        filterGtf
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { WorkflowRiboseq.checkMaxContigSize(it, log) }
    }

    //
    // Create input channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            WorkflowRiboseq.validateInput(it)
        }
        .set { ch_fastq }

    //
    // SUBWORKFLOW: preprocess reads for RNA-seq. Includes trimming,
    // contaminant removal, strandedness inference
    //

    PREPROCESS_RNASEQ (
        ch_fastq,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !params.salmon_index && !('salmon' in prepareToolIndices),
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read
    )
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_RNASEQ.out.multiqc_files)
    ch_versions      = ch_versions.mix(PREPROCESS_RNASEQ.out.versions)

    //
    // SUBWORKFLOW: align with STAR, produce both genomic and transcriptomic
    // alignments and run BAM_SORT_STATS_SAMTOOLS for each
    //

    FASTQ_ALIGN_STAR(
        PREPROCESS_RNASEQ.out.reads,
        PREPARE_GENOME.out.star_index.map { [ [:], it ] },
        PREPARE_GENOME.out.gtf.map { [ [:], it ] },
        params.star_ignore_sjdbgtf,
        '',
        params.seq_center ?: '',
        PREPARE_GENOME.out.fasta.map { [ [:], it ] },
        PREPARE_GENOME.out.transcript_fasta.map { [ [:], it ] }
    )

    ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
    ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
    ch_transcriptome_bam       = FASTQ_ALIGN_STAR.out.bam_transcript
    ch_transcriptome_bam_index = FASTQ_ALIGN_STAR.out.bai_transcript
    ch_versions                = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

    ch_multiqc_files = ch_multiqc_files
        .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //

    if (params.with_umi) {

        // Deduplicate genome BAM file before downstream analysis

        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
            ch_genome_bam.join(ch_genome_bai, by: [0]),
            params.umitools_dedup_stats
        )
        ch_genome_bam       = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
        ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai

        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats.collect{it[1]})
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat.collect{it[1]})
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats.collect{it[1]})

        // Deduplicate transcriptome BAM file before downstream analysis

        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
            ch_transcriptome_bam.join(ch_transcriptome_bai, by: [0]),
            params.umitools_dedup_stats
        )
        ch_transcriptome_bam       = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
        ch_transcriptome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bai

        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.stats.collect{it[1]})
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.flagstat.collect{it[1]})
            .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.idxstats.collect{it[1]})

        // Prepare trancriptome BAM for Salmon. This requires a Samtools name
        // sort and a specific umittools command

        SAMTOOLS_SORT (
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
        )

        // Only run prepare_for_rsem.py on paired-end BAM files
        SAMTOOLS_SORT.out.bam
            .branch { meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
            }
            .set { ch_umitools_dedup_bam }

        // Fix paired-end reads in name sorted BAM file
        // See: https://github.com/nf-core/rnaseq/issues/828
        UMITOOLS_PREPAREFORSALMON (
            ch_umitools_dedup_bam.paired_end.map{meta, bam -> [meta, bam, []]}
        )
        ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

        ch_umitools_dedup_bam.single_end
            .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
            .set { ch_transcriptome_bam_for_salmon }
    }


    //
    // Compile software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRiboseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRiboseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()

    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
