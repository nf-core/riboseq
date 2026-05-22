/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS                                                 } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main'
include { FASTQ_EQUALISE_READ_LENGTHS                                                          } from '../../subworkflows/local/fastq_equalise_read_lengths'
include { BAM_DEDUP_UMI      } from '../../subworkflows/nf-core/bam_dedup_umi'
include { FASTQ_ALIGN_STAR   } from '../../subworkflows/nf-core/fastq_align_star'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                                              } from '../../modules/nf-core/multiqc/main'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../modules/nf-core/umitools/prepareforrsem'
include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_RIBOSEQ         } from '../../modules/nf-core/ribotish/quality'
include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_TISEQ           } from '../../modules/nf-core/ribotish/quality'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_INDIVIDUAL      } from '../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_ALL             } from '../../modules/nf-core/ribotish/predict'
include { RIBOTRICER_PREPAREORFS                               } from '../../modules/nf-core/ribotricer/prepareorfs'
include { RIBOTRICER_DETECTORFS                                } from '../../modules/nf-core/ribotricer/detectorfs'
include { RIBOCODE_GTFUPDATE                                   } from '../../modules/nf-core/ribocode/gtfupdate'
include { RIBOCODE_PREPARE                                     } from '../../modules/nf-core/ribocode/prepare'
include { RIBOCODE_METAPLOTS                                   } from '../../modules/nf-core/ribocode/metaplots'
include { RIBOCODE_RIBOCODE                                    } from '../../modules/nf-core/ribocode/ribocode'
include { ANOTA2SEQ_ANOTA2SEQRUN                               } from '../../modules/nf-core/anota2seq/anota2seqrun'
include { DESEQ2_DELTATE                                       } from '../../modules/local/deseq2/deltate'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_STAR_SALMON    } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_PSEUDO_TE      } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { RIBOWALTZ                                            } from '../../modules/nf-core/ribowaltz/main'
include { PLASTID_METAGENE_GENERATE                            } from '../../modules/nf-core/plastid/metagene_generate/main'
include { PLASTID_PSITE                                        } from '../../modules/nf-core/plastid/psite/main'
include { PLASTID_MAKE_WIGGLE                                  } from '../../modules/nf-core/plastid/make_wiggle/main'
include { QUANTIFY_INFRAME_PSITE_PLASTID                       } from '../../modules/local/quantify_inframe_psite_plastid'
include { GAWK as GTF_TO_INFRAME_PSITES                        } from '../../modules/nf-core/gawk'
include { GAWK as REPLACE_RIBOSEQ_COUNTS_IN_MATRIX             } from '../../modules/nf-core/gawk'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SPLIT_BY_STRAND       } from '../../modules/nf-core/samtools/view'
include { BEDTOOLS_GENOMECOV                                   } from '../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG                                } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap         } from 'plugin/nf-schema'
include { samplesheetToList        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText   } from '../../subworkflows/local/utils_nfcore_riboseq_pipeline'
include { validateInputSamplesheet } from '../../subworkflows/local/utils_nfcore_riboseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// A filter for samtools view which splits alignments by first-of-pair strand,
// taking into consideration the strandedness of the library. Used by SAMTOOLS_VIEW_SPLIT_BY_STRAND.
def getStrandFilter(strandedness, strand) {
    def sameOrientation = (strand == 'forward') == (strandedness == 'forward')
    sameOrientation
        ? "-e '((flag.read1 || !flag.paired) && !flag.reverse) || (flag.read2 &&  flag.reverse)'"
        : "-e '((flag.read1 || !flag.paired) &&  flag.reverse) || (flag.read2 && !flag.reverse)'"
}

workflow RIBOSEQ {

    take:
    ch_samplesheet      // channel: path(sample_sheet.csv)
    ch_contrasts_file   // channel: path(contrasts.csv)
    ch_versions         // channel: [ path(versions.yml) ]
    ch_fasta            // channel: path(genome.fasta)
    ch_gtf              // channel: path(genome.gtf)              - full multi-isoform, used for genome-guided alignment
    ch_canonical_gtf    // channel: path(canonical.gtf)           - one-transcript-per-gene backbone for ORF calling, P-site calibration, DTE
    ch_fai              // channel: path(genome.fai)
    ch_chrom_sizes      // channel: path(genome.sizes)
    ch_transcript_fasta // channel: path(transcript.fasta)
    ch_star_index       // channel: path(star/index/)
    ch_salmon_index     // channel: path(salmon/index/)
    ch_salmon_index_te  // channel: path(salmon_te/index/) - for TE pseudo-alignment
    ch_bbsplit_index    // channel: path(bbsplit/index/)
    ch_rrna_fastas      // channel: path(fasta)
    ch_sortmerna_index  // channel: path(sortmerna/index/)
    ch_bowtie2_index    // channel: path(bowtie2/index/) for rRNA removal

    main:

    //
    // Collect versions from topic channel (for modules that emit versions via topics)
    //
    def topic_versions = channel.topic('versions')
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Check rRNA databases for sortmerna/bowtie2 (ribodetector doesn't need a database)
    if (params.remove_ribo_rna && params.ribo_removal_tool != 'ribodetector') {
        ch_ribo_db = file(params.ribo_database_manifest)
        if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
    } else {
        ch_ribo_db = Channel.empty()
    }

    // Check if file with list of fastas is provided when running BBSplit
    if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
        ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list)
        if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
    }

    // Check alignment parameters
    def prepareToolIndices  = []
    if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
    if (params.remove_ribo_rna) { prepareToolIndices << 'sortmerna' }
    if (!params.skip_alignment) { prepareToolIndices << params.aligner }

    ch_multiqc_files = Channel.empty()

    //
    // Create input channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
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
            validateInputSamplesheet(it)
        }
        .set { ch_fastq }

    //
    // SUBWORKFLOW: preprocess reads for RNA-seq. Includes trimming,
    // contaminant removal, strandedness inference
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    // Determine if we need to build rRNA removal indexes
    def make_sortmerna_index = !params.sortmerna_index && params.remove_ribo_rna && params.ribo_removal_tool == 'sortmerna'
    def make_bowtie2_index   = params.remove_ribo_rna && params.ribo_removal_tool == 'bowtie2'

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,                                   // ch_reads
        ch_fasta,                                   // ch_fasta
        ch_transcript_fasta,                        // ch_transcript_fasta
        ch_gtf,                                     // ch_gtf
        ch_salmon_index,                            // ch_salmon_index
        ch_sortmerna_index,                         // ch_sortmerna_index
        ch_bowtie2_index,                           // ch_bowtie2_index
        ch_bbsplit_index,                           // ch_bbsplit_index
        ch_rrna_fastas,                             // ch_rrna_fastas
        params.skip_bbsplit,                        // skip_bbsplit
        params.skip_fastqc || params.skip_qc,       // skip_fastqc
        params.skip_trimming,                       // skip_trimming
        params.skip_umi_extract,                    // skip_umi_extract
        params.skip_linting,                        // skip_linting
        !salmon_index_available,                    // make_salmon_index
        make_sortmerna_index,                       // make_sortmerna_index
        make_bowtie2_index,                         // make_bowtie2_index
        params.trimmer,                             // trimmer
        params.min_trimmed_reads,                   // min_trimmed_reads
        params.save_trimmed,                        // save_trimmed
        false,                                      // fastp_merge
        params.remove_ribo_rna,                     // remove_ribo_rna
        params.ribo_removal_tool,                   // ribo_removal_tool
        params.with_umi,                            // with_umi
        params.umi_discard_read,                    // umi_discard_read
        false,                                      // save_merged_fastq
        params.stranded_threshold,                  // stranded_threshold
        params.unstranded_threshold                 // unstranded_threshold
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)

    //
    // SUBWORKFLOW: Equalise RNA-seq read lengths to match Ribo-seq read lengths
    //

    // Normalize samplesheet-derived meta fields (convert empty lists to null)
    // Only modify fields that exist in the meta
    ch_reads_preprocessed = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
        .map { meta, reads ->
            def updates = [:]
            if (meta.containsKey('trim_length')) {
                updates.trim_length = meta.trim_length instanceof List ? (meta.trim_length ? meta.trim_length[0] : null) : meta.trim_length
            }
            if (meta.containsKey('pair')) {
                updates.pair = meta.pair instanceof List ? (meta.pair ? meta.pair[0] : null) : meta.pair
            }
            return [ meta + updates, reads ]
        }

    ch_reads_for_alignment = ch_reads_preprocessed

    if (params.equalise_read_lengths) {
        FASTQ_EQUALISE_READ_LENGTHS(
            ch_reads_preprocessed,
            params.equalise_read_lengths_target
        )
        ch_reads_for_alignment = FASTQ_EQUALISE_READ_LENGTHS.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EQUALISE_READ_LENGTHS.out.riboseq_stats.collect{it[1]})
    }

    //
    // SUBWORKFLOW: align with STAR, produce both genomic and transcriptomic
    // alignments and run BAM_SORT_STATS_SAMTOOLS for each
    //

    FASTQ_ALIGN_STAR(
        ch_reads_for_alignment,
        ch_star_index.map { [ [:], it ] },
        ch_gtf.map { [ [:], it ] },
        params.star_ignore_sjdbgtf,
        ch_fasta.map { [ [:], it ] },
        ch_transcript_fasta.map { [ [:], it ] }
    )

    ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
    ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
    ch_transcriptome_bam       = FASTQ_ALIGN_STAR.out.orig_bam_transcript
    ch_transcriptome_bai       = FASTQ_ALIGN_STAR.out.bai_transcript

    ch_multiqc_files = ch_multiqc_files
        .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //

    if (params.with_umi) {

        BAM_DEDUP_UMI(
            ch_genome_bam.join(ch_genome_bam_index, by: [0]),
            ch_fasta.map { [ [:], it ] },
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            params.bam_csi_index,
            ch_transcriptome_bam,
            ch_transcript_fasta.map { [ [:], it ] },
            params.umitools_dedup_primary_only
        )

        ch_genome_bam        = BAM_DEDUP_UMI.out.bam
        ch_transcriptome_bam = BAM_DEDUP_UMI.out.transcriptome_bam
        ch_genome_bam_index  = BAM_DEDUP_UMI.out.bai

        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_DEDUP_UMI.out.multiqc_files)
    }

    //
    // Generate coverage tracks
    //

    if (!params.skip_coverage_tracks) {

        // When protocol is stranded, split BAMs by mate1 strand
        ch_split_by_strand = ch_genome_bam
            .join(ch_genome_bam_index, by: [0])
            .filter { meta, bam, bai -> meta.strandedness in ['forward', 'reverse'] }
        SAMTOOLS_VIEW_SPLIT_BY_STRAND(
            ch_split_by_strand
                .flatMap { meta, bam, bai ->
                    ['forward', 'reverse'].collect { strand ->
                        def strand_filter = getStrandFilter(meta.strandedness, strand)
                        [meta + [strand: strand, strand_filter: strand_filter], bam, bai]
                    }
                },
            [[], [], []],  // No reference fasta/fai
            [],            // No qname file
            []             // No index format
        )

        // Create bedgraph tracks
        BEDTOOLS_GENOMECOV(
            SAMTOOLS_VIEW_SPLIT_BY_STRAND.out.bam
                .map { meta, bam -> [meta, bam, 1] }
                .mix(ch_genome_bam
                    .filter { meta, bam -> meta.strandedness == 'unstranded' }
                    .map { meta, bam -> [meta + [strand: 'unstranded'], bam, 1] }
                ),
            ch_fai,
            'bedgraph',
            false
        )
        // Convert bedgraphs to bigWig
        UCSC_BEDGRAPHTOBIGWIG(
            BEDTOOLS_GENOMECOV.out.genomecov,
            ch_fai
        )

    }

    //
    // Take the riboseq samples and route to ribotish
    //

    ch_genome_bam
        .branch { meta, bam ->
            riboseq: meta.sample_type == 'riboseq'
                return [ meta, bam ]
            tiseq: meta.sample_type == 'tiseq'
                return [ meta, bam ]
            rnaseq: meta.sample_type == 'rnaseq'
                return [ meta, bam ]
        }
        .set{
            ch_genome_bam_by_type
        }

    ch_bams_for_analysis = ch_genome_bam_by_type.riboseq.join(ch_genome_bam_index)
    // Use the canonical (one-transcript-per-gene) annotation backbone for ORF
    // calling, P-site calibration and DTE; the full `ch_gtf` is reserved for
    // genome-guided alignment.
    ch_fasta_gtf = ch_fasta.combine(ch_canonical_gtf).map{ fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }.first()

    if (!params.skip_ribotish){
        RIBOTISH_QUALITY_RIBOSEQ(
            ch_bams_for_analysis,
            ch_canonical_gtf.map { [ [:], it ] }.first()
        )
        ch_versions      = ch_versions.mix(RIBOTISH_QUALITY_RIBOSEQ.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RIBOTISH_QUALITY_RIBOSEQ.out.distribution.collect{it[1]})

        ribotish_predict_inputs = ch_bams_for_analysis
            .join(RIBOTISH_QUALITY_RIBOSEQ.out.offset)
            .multiMap{ meta, bam, bai, offset ->
                bam: [ meta, bam, bai ]
                offset: [ meta, offset ]
            }

        RIBOTISH_PREDICT_INDIVIDUAL(
            ribotish_predict_inputs.bam,
            [[:],[],[]],
            ch_fasta_gtf,
            [[:],[]],
            ribotish_predict_inputs.offset,
            [[:],[]]
        )
        ch_versions = ch_versions.mix(RIBOTISH_PREDICT_INDIVIDUAL.out.versions)

        RIBOTISH_PREDICT_ALL(
            ribotish_predict_inputs.bam.map{meta, bam, bai -> [[id:'allsamples'], bam, bai]}.groupTuple(),
            [[:],[],[]],
            ch_fasta_gtf,
            [[:],[]],
            ribotish_predict_inputs.offset.map{meta, offset -> [[id:'allsamples'], offset]}.groupTuple(),
            [[:],[]]
        )
        ch_versions = ch_versions.mix(RIBOTISH_PREDICT_ALL.out.versions)
    }

    if (params.run_ribotricer){
        log.warn "Ribotricer is enabled via --run_ribotricer. Benchmark data (FK/NGB, May 2026) found its ORF-score column is rank-unstable across biological replicates (mean Spearman 0.288 vs Jaccard 0.770). Its binary calls are usable, but do not rely on its scores as the primary ranking source; the cross-caller rank aggregation will exclude them."

        RIBOTRICER_PREPAREORFS(
            ch_fasta_gtf
        )
        ch_versions = ch_versions.mix(RIBOTRICER_PREPAREORFS.out.versions)

        RIBOTRICER_DETECTORFS(
            ch_bams_for_analysis,
            RIBOTRICER_PREPAREORFS.out.candidate_orfs
        )
        ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS.out.versions)
    }

    //
    // Dynamic ORF-caller set for cross-caller agreement (issue #07).
    // The enabled list reflects which callers ran at runtime; the agreement
    // threshold and rank-aggregation set are derived from it so the logic
    // works whether 2 (default) or 3 callers are active.
    //
    def enabled_orf_callers = []
    if (!params.skip_ribotish)   { enabled_orf_callers << 'ribotish' }
    if (!params.skip_ribocode)   { enabled_orf_callers << 'ribocode' }
    if ( params.run_ribotricer)  { enabled_orf_callers << 'ribotricer' }

    // Ribotricer contributes binary calls only; its scores are excluded from
    // the cross-caller rank aggregation due to known rank instability.
    def rank_aggregation_callers  = enabled_orf_callers - 'ribotricer'
    // Strict-majority of enabled callers (floor(N/2)+1): N=2 -> 2 (both must
    // agree), N=3 -> 2 (majority). Adapts as the caller set grows.
    def orf_agreement_min_callers = enabled_orf_callers
        ? enabled_orf_callers.size().intdiv(2) + 1
        : 0
    ch_enabled_orf_callers      = Channel.value(enabled_orf_callers)
    ch_rank_aggregation_callers = Channel.value(rank_aggregation_callers)

    if (!params.skip_ribocode){
        // RiboCode requires transcriptome BAMs
        ch_transcriptome_bams_for_ribocode = ch_transcriptome_bam
            .branch { meta, bam ->
                riboseq: meta.sample_type == 'riboseq'
                    return [ meta, bam ]
            }
            .riboseq

        // Step 1: Update GTF annotation
        RIBOCODE_GTFUPDATE(
            ch_gtf.map { [ [:], it ] }.first()
        )

        // Step 2: Prepare annotation files
        RIBOCODE_PREPARE(
            ch_fasta.map { [ [:], it ] }.first(),
            RIBOCODE_GTFUPDATE.out.gtf
        )

        // Step 3: Generate metaplots and config for each sample
        RIBOCODE_METAPLOTS(
            ch_transcriptome_bams_for_ribocode,
            RIBOCODE_PREPARE.out.annotation
        )

        // Step 4: Run RiboCode ORF detection
        // Join BAMs with their corresponding config files by meta
        ch_ribocode_inputs = ch_transcriptome_bams_for_ribocode
            .join(RIBOCODE_METAPLOTS.out.config)

        RIBOCODE_RIBOCODE(
            ch_ribocode_inputs.map { meta, bam, config -> [ meta, bam ] },
            RIBOCODE_PREPARE.out.annotation,
            ch_ribocode_inputs.map { meta, bam, config -> [ meta, config ] }
        )
    }


    //
    // Get P-sites and P-site diagnostics with riboWaltz
    //

    // Keep only riboseq transcriptome BAMs for riboWaltz
    ch_transcriptome_bam
        .branch { meta, bam ->
            riboseq: meta.sample_type == 'riboseq'
                return [ meta, bam ]
            tiseq: meta.sample_type == 'tiseq'
                return [ meta, bam ]
            rnaseq: meta.sample_type == 'rnaseq'
                return [ meta, bam ]
        }
        .set { ch_transcriptome_bam_by_type }

    if (!params.skip_ribowaltz) {
        RIBOWALTZ(
            ch_transcriptome_bam_by_type.riboseq,
            ch_gtf.map { [ [:], it ] },
            ch_fasta.map { [ [:], it ] })

        ch_versions = ch_versions.mix(RIBOWALTZ.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RIBOWALTZ.out.ribowaltz_qc_data.collect{it[1]}.ifEmpty([]))
    }

    if (!params.skip_plastid) {

        PLASTID_METAGENE_GENERATE(ch_canonical_gtf.map { [ [:], it ] })
        ch_versions = ch_versions.mix(PLASTID_METAGENE_GENERATE.out.versions)

        PLASTID_PSITE(
            ch_bams_for_analysis,
            PLASTID_METAGENE_GENERATE.out.rois_txt
        )
        ch_versions = ch_versions.mix(PLASTID_PSITE.out.versions)

        PLASTID_MAKE_WIGGLE(
            ch_bams_for_analysis.join(PLASTID_PSITE.out.p_offsets, by: [0]),
            "fiveprime_variable"
        )
        ch_versions = ch_versions.mix(PLASTID_MAKE_WIGGLE.out.versions)

    }

    //
    // SUBWORKFLOW: Count reads from BAM alignments using Salmon
    //

    // Salmon transcriptome quantification uses the full GTF: tx2gene must match
    // the transcript fasta the Salmon index was built against.
    QUANTIFY_STAR_SALMON (
        ch_samplesheet.map { [ [:], it ] },
        ch_transcriptome_bam,
        [],
        ch_transcript_fasta,
        ch_gtf,
        params.gtf_group_features,
        params.gtf_extra_attributes,
        'salmon',
        true,
        params.salmon_quant_libtype ?: '',
        null,
        null
    )
    ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_STAR_SALMON.out.multiqc.collect{it[1]}.ifEmpty([]))

    //
    // SUBWORKFLOW: Pseudo-alignment quantification for TE analysis (when enabled)
    // Uses direct Salmon pseudo-alignment with a lower k-mer index optimized for short Ribo-seq reads
    //

    ch_te_counts = QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled  // Default: use alignment-based counts

    if (params.te_quantification_method == 'pseudo' && ch_contrasts_file) {
        // Filter reads to only riboseq and rnaseq for TE pseudo-alignment
        ch_reads_for_te = ch_reads_for_alignment
            .filter { meta, reads -> meta.sample_type in ['riboseq', 'rnaseq'] }

        QUANTIFY_PSEUDO_TE (
            ch_samplesheet.map { [ [:], it ] },
            ch_reads_for_te,
            ch_salmon_index_te,
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            false,  // alignment_mode = false (pseudo-alignment from reads)
            params.salmon_quant_libtype ?: '',
            null,
            null
        )
        ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_TE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_TE.out.multiqc.collect{it[1]}.ifEmpty([]))
        ch_te_counts = QUANTIFY_PSEUDO_TE.out.counts_gene_length_scaled
    }

    if (params.te_quantification_method == 'plastid_psite' && !params.skip_plastid) {
        // Convert GTF CDS segments to in-frame p-site positions
        GTF_TO_INFRAME_PSITES(
            ch_canonical_gtf.map { gtf -> [ [id: gtf.baseName, feature: 'gene'], gtf ] },
            file("${projectDir}/bin/gtf_to_inframe_psites.awk"),
            false
        )

        ch_inframe_psites = GTF_TO_INFRAME_PSITES.out.output.first()

        // Run p-site quantification per sample
        ch_psite_tracks = PLASTID_MAKE_WIGGLE.out.tracks
            .map { meta, tracks -> [meta, tracks[0], tracks[1]] }

        QUANTIFY_INFRAME_PSITE_PLASTID(
            ch_psite_tracks,
            ch_inframe_psites
        )

        // Merge per-sample p-site counts into a single file
        ch_psite_counts_merged = QUANTIFY_INFRAME_PSITE_PLASTID.out.counts
            .collectFile(name: 'gene_inframe_psite_counts.tsv') { meta, file -> file }
            .map { file -> [ [:], file ] }

        REPLACE_RIBOSEQ_COUNTS_IN_MATRIX(
            ch_psite_counts_merged
                .combine(QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map{ meta, counts -> counts })
                .map { meta, psite_counts, salmon_counts -> [meta, [psite_counts, salmon_counts]] },
            file("${projectDir}/bin/replace_riboseq_counts_in_matrix.awk"),
            false
        )
        ch_te_counts = REPLACE_RIBOSEQ_COUNTS_IN_MATRIX.out.output
        ch_versions = ch_versions.mix(QUANTIFY_INFRAME_PSITE_PLASTID.out.versions)
    }

    //
    // Do a translational efficiency analysis where contrasts are supplied
    //

    if (ch_contrasts_file){

        ch_contrasts = ch_contrasts_file
            .splitCsv ( header:true, sep:',' )
            .map{[it, it.variable, it.reference, it.target]}

        ch_samplesheet_matrix = ch_te_counts
            .combine(ch_samplesheet)
            .map{[it[0], it[2], it[1]]}
            .first()

        if (params.translational_efficiency_method == 'anota2seq') {
            ANOTA2SEQ_ANOTA2SEQRUN(
                ch_contrasts,
                ch_samplesheet_matrix
            )
            ch_versions = ch_versions.mix(ANOTA2SEQ_ANOTA2SEQRUN.out.versions)
        }

        if (params.translational_efficiency_method == 'deltate') {
            DESEQ2_DELTATE(
                ch_contrasts,
                ch_samplesheet_matrix
            )
            ch_versions = ch_versions.mix(DESEQ2_DELTATE.out.versions)
        }
    }

    //
    // Collate and save software versions
    // Combines traditional versions.yml files with versions emitted via topic channels
    //
    ch_versions = ch_versions.filter{it != null}

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_riboseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        ch_name_replacements = ch_fastq
            .map{ meta, reads ->
                def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
                def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
                if (reads[0][1] ){
                    def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
                    def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
                    def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
                    return [ name1, name2, fastqcnames1, fastqcnames2 ]
                } else{
                    return [ name1, fastqcnames ]
                }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC(
            ch_multiqc_files.flatten().collect()
                .combine(ch_name_replacements)
                .map { args ->
                    def replace_names = args[-1]
                    def files = args[0..-2]
                    def mqc_config = params.multiqc_config
                        ? file(params.multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
                    [
                        [id: 'riboseq'],
                        files,
                        mqc_config,
                        params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                        replace_names,
                        [],
                    ]
                }
        )
    ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    } else {
        ch_multiqc_report = Channel.empty()
    }

    emit:
    multiqc_report = ch_multiqc_report   // channel: /path/to/multiqc_report.html
    versions       = ch_versions         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
