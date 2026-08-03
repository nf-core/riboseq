//
// Subworkflow with functionality specific to the nf-core/riboseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/riboseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/riboseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet, params.with_umi)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    removedParamsError()
    dotseqPrerequisitesError()
    kallistoPrerequisitesError()

    // --extended_orf_analysis is a no-op without a novel-transcript source
    // (StringTie or a user-supplied --novel_gtf); warn rather than error so
    // flags can be composed incrementally.
    def novel_source_configured = !params.skip_stringtie || params.novel_gtf
    if (params.extended_orf_analysis && !novel_source_configured) {
        log.warn "--extended_orf_analysis is enabled but no novel-transcript source is configured (--skip_stringtie is true and --novel_gtf is unset). The flag has no effect; ORF callers will run against the canonical GTF as usual."
    }
    if (params.extended_orf_analysis && novel_source_configured && params.skip_plastid) {
        log.warn "--extended_orf_analysis is enabled but --skip_plastid is true. ORF-level P-site quantification needs the plastid wiggle tracks and will be skipped; the ORF catalogue will still be built."
    }
}

//
// Exit pipeline on removed parameters. Schema validation only warns about
// unrecognised parameters, too quiet for flags whose purpose was to suppress
// work the pipeline does regardless.
//
def removedParamsError() {
    def removed = [
        'min_mapped_reads'     : 'it was never applied to any sample',
        'skip_pseudo_alignment': 'pseudo-alignment is selected with --te_quantification_method and --pseudo_aligner',
        'skip_alignment'       : 'the pipeline has no index-only mode; every downstream stage needs the alignment',
    ]

    def supplied = removed.findAll { name, _reason -> params.containsKey(name) }
    if (supplied) {
        error("The following parameters have been removed:\n" + supplied.collect { name, reason -> "  --${name}: ${reason}" }.join('\n'))
    }
}

//
// Exit pipeline if kallisto is selected as the pseudo-aligner without what it
// needs: fragment length stats (single-end libraries) and a valid k-mer size.
//
def kallistoPrerequisitesError() {
    if (params.pseudo_aligner != 'kallisto' || params.te_quantification_method != 'pseudo') return

    if (!params.kallisto_quant_fraglen || !params.kallisto_quant_fraglen_sd) {
        error("--pseudo_aligner kallisto requires --kallisto_quant_fraglen and --kallisto_quant_fraglen_sd, which kallisto needs to quantify single-end libraries. Set both, or use --pseudo_aligner salmon.")
    }

    if (!params.kallisto_index && (params.pseudo_aligner_kmer_size % 2 == 0 || params.pseudo_aligner_kmer_size > 31)) {
        error("--pseudo_aligner_kmer_size ${params.pseudo_aligner_kmer_size} is invalid for kallisto index building: kallisto requires an odd k-mer size no greater than 31. Set --pseudo_aligner_kmer_size to an odd value <= 31, or supply a pre-built --kallisto_index.")
    }
}

//
// Exit pipeline if --translational_efficiency_method dotseq is selected
// without the ORF-level prerequisites that DOTSeq needs.
//
def dotseqPrerequisitesError() {
    if (!('dotseq' in params.translational_efficiency_method.tokenize(',')*.trim())) return

    def novel_source_configured = !params.skip_stringtie || params.novel_gtf
    def extended_orf_active     = params.extended_orf_analysis && novel_source_configured
    def any_caller_enabled      = !params.skip_ribotish || !params.skip_ribocode || params.run_ribotricer || params.run_rpbp || params.run_price

    if (!extended_orf_active || !any_caller_enabled || params.skip_plastid || !params.contrasts) {
        def missing = []
        if (!params.extended_orf_analysis) missing << "--extended_orf_analysis true"
        else if (!novel_source_configured) missing << "a novel-transcript source (set --novel_gtf or leave --skip_stringtie false)"
        if (!any_caller_enabled)           missing << "at least one ORF caller (do not skip both ribocode and ribotish, or opt into ribotricer / rpbp / price)"
        if (params.skip_plastid)           missing << "--skip_plastid false"
        if (!params.contrasts)             missing << "--contrasts"
        error("--translational_efficiency_method dotseq runs only at the ORF level and requires: ${missing.join('; ')}. Pick anota2seq or deltate for the gene-level path, or wire up the ORF prerequisites.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input, default_with_umi) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    def umi_status_ok = metas.collect { meta -> resolveSampleUmi(meta, default_with_umi) }.unique().size == 1
    if (!umi_status_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same with_umi value: ${metas[0].id}")
    }

    def clean_meta = new LinkedHashMap(metas[0])
    clean_meta.remove('with_umi')
    return [ clean_meta, fastqs ]
}

def resolveSampleUmi(meta, default_with_umi) {
    def sample_with_umi = meta.with_umi
    if (sample_with_umi instanceof List) {
        sample_with_umi = sample_with_umi ? sample_with_umi[0] : null
    }
    if (sample_with_umi == null || sample_with_umi == '') {
        sample_with_umi = default_with_umi
    }
    if (!(sample_with_umi instanceof Boolean)) {
        error("Please check input samplesheet -> with_umi must be true or false for sample: ${meta.id}")
    }
    return sample_with_umi
}

def samplesheetUmiSampleIds(samplesheet_rows, default_with_umi) {
    samplesheet_rows
        .groupBy { meta, _fastq_1, _fastq_2 -> meta.id }
        .findAll { sample_id, rows ->
            def statuses = rows.collect { meta, _fastq_1, _fastq_2 -> resolveSampleUmi(meta, default_with_umi) }.unique()
            if (statuses.size() != 1) {
                error("Please check input samplesheet -> Multiple runs of a sample must have the same with_umi value: ${sample_id}")
            }
            statuses[0]
        }
        .keySet() as Set
}

def trimFailuresMultiqcTsv(trim_read_counts, min_trimmed_reads) {
    def failures = trim_read_counts.findAll { meta, num_reads -> num_reads <= min_trimmed_reads.toFloat() }
        .collect { meta, num_reads -> "${meta.id}\t${num_reads}" }

    failures ? "Sample\tReads after trimming\n${failures.join('\n')}" : ''
}

def samplesheetNeedsSalmonForStrandedness(input) {
    samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .any { meta, _fastq_1, _fastq_2 -> meta.strandedness == 'auto' }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function to generate an error if contigs in genome fasta file > 512 Mbp
//
def checkMaxContigSize(fai_file) {
    def max_size = 512000000
    fai_file.eachLine { line ->
        def lspl  = line.split('\t')
        def chrom = lspl[0]
        def size  = lspl[1]
        if (size.toInteger() > max_size) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
                "  ${chrom}: ${size}\n\n" +
                "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
                "  Please see:\n" +
                "  https://github.com/nf-core/rnaseq/issues/744\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            error(error_string)
        }
    }
}
