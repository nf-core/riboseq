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

WorkflowRiboseq.initialise(params, log)

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                                                        } from '../modules/nf-core/fastqc/main'
include { HISAT2_EXTRACTSPLICESITES                                                     } from '../modules/nf-core/hisat2/extractsplicesites/main' 
include { HISAT2_BUILD as HISAT2_BUILD_rRNA; HISAT2_BUILD as HISAT2_BUILD_transcriptome } from '../modules/nf-core/hisat2/build/main'                               
include { RSEM_PREPAREREFERENCE                                                         } from '../modules/nf-core/rsem/preparereference/main'                                                  
include { UMITOOLS_EXTRACT                                                              } from '../modules/nf-core/umitools/extract/main'                                                            
include { CUTADAPT                                                                      } from '../modules/nf-core/cutadapt/main'
include { HISAT2_ALIGN as HISAT2_ALIGN_rRNA; HISAT2_ALIGN as HISAT2_ALIGN_transcriptome } from '../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_SORT                                                                 } from '../modules/nf-core/samtools/sort/main' 
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_transcriptome; SAMTOOLS_INDEX as SAMTOOLS_INDEX_umi_dedup } from '../modules/nf-core/samtools/index/main'                                                                
include { UMITOOLS_DEDUP                                                                } from '../modules/nf-core/umitools/dedup/main'
include { BEDTOOLS_BAMTOBED                                                             } from '../modules/nf-core/bedtools/bamtobed/main'                          
include { BEDTOOLS_GENOMECOV                                                            } from '../modules/nf-core/bedtools/genomecov/main'                                  
include { SUBREAD_FEATURECOUNTS                                                         } from '../modules/nf-core/subread/featurecounts/main'                                                  
include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RIBOSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    
    ch_fastq = Channel.fromSamplesheet("input")
    ch_fastq.view()

    // Creating channels for HISAT2_BUILD and HISAT2_ALIGN
    ch_transcriptome_fasta = Channel.fromPath(params.transcriptome_fasta)
    ch_genome_fasta = Channel.fromPath(params.genome_fasta)
    ch_gtf = Channel.fromPath(params.gtf)
    ch_rRNA_fasta = Channel.fromPath(params.rRNA_fasta)
    //ch_gtf_modified = Channel.fromPath(params.gtf_modified)
    //ch_transcriptome_modified = Channel.fromPath(params.transcriptome_modified)
    

    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    // MODULE: Run Hisat2_extractsplicesites
    //
    ch_splicesites = HISAT2_EXTRACTSPLICESITES (
        ch_gtf.map { [ [:], it ] } 
        )



    //MODULE: Run Hisat2_Build for rRNA
    //
    HISAT2_BUILD_rRNA (
        ch_rRNA_fasta.map { [ [:], it ] },
        [[],[]],
        [[],[]]
    )
    

    // MODULE: Run Hisat2_Build for transcriptome
    //
    HISAT2_BUILD_transcriptome (
        ch_transcriptome_fasta.map { [ [:], it ] },
        ch_gtf.map { [ [:], it ] },
        ch_splicesites.txt
    )



    // MODULE: Run cutadapt
    //
    CUTADAPT (
        ch_fastq
    )

    // MODULE: Run UMITOOLS_EXTRACT
    //
    with_umi = params.with_umi
    skip_umi_extract = params.skip_umi_extract
    if (with_umi && !skip_umi_extract) {
    UMITOOLS_EXTRACT (
       CUTADAPT.out.reads
    )
    }


    // MODULE: Run Hisat2_Align for rRNA 
    //

    if (!params.skip_alignment && params.aligner == 'hisat2') {
    HISAT2_ALIGN_rRNA (
        UMITOOLS_EXTRACT.out.reads,
        HISAT2_BUILD_rRNA.out.index,
        [[],[]]
    )
    }
 
    // MODULE: Run Hisat2_Align for transcriptome
    //
    HISAT2_ALIGN_transcriptome (
        HISAT2_ALIGN_rRNA.out.fastq,
        HISAT2_BUILD_transcriptome.out.index,
        ch_splicesites.txt
    )

    // MODULE: Run SAMTOOLS_SORT
    //
    SAMTOOLS_SORT (
        HISAT2_ALIGN_transcriptome.out.bam
    )

    // MODULE: Run SAMTOOLS_INDEX
    //
    SAMTOOLS_INDEX_transcriptome (
        SAMTOOLS_SORT.out.bam
    )

    // MODULE: Run UMITOOLS_DEDUP 
    // 
    
    ch_transcriptome_sorted_bam = SAMTOOLS_SORT.out.bam
    ch_transcriptome_sorted_bai = SAMTOOLS_INDEX_transcriptome.out.bai

    // Deduplicate genome BAM file before downstream analysis
    if (params.with_umi) {
    UMITOOLS_DEDUP (
        ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
            params.umitools_dedup_stats
    )
    }
    
    // MODULE: Run SAMTOOLS_INDEX
    //
    SAMTOOLS_INDEX_umi_dedup (
        UMITOOLS_DEDUP.out.bam
    )

    // MODULE: Run BEDTOOLS_BAMTOBED
    //
    BEDTOOLS_BAMTOBED (
        UMITOOLS_DEDUP.out.bam
    )


    // MODULE: Run BEDTOOLS_GENOMECOV
    //
    //BEDTOOLS_GENOMECOV (
       // UMITOOLS_DEDUP.out.bam
    //)

    // MODULE: Run Featurecoun
    //
    SUBREAD_FEATURECOUNTS (
        UMITOOLS_DEDUP.out.bam
            .combine(ch_gtf)
    )


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRiboseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRiboseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
