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
include { BAM_DEDUP_UMI       } from '../../subworkflows/nf-core/bam_dedup_umi'
include { FASTQ_ALIGN_STAR    } from '../../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_STAR as FASTQ_ALIGN_STAR_HYBRID } from '../../subworkflows/nf-core/fastq_align_star'
include { BAM_STRINGTIE_MERGE      } from '../../subworkflows/nf-core/bam_stringtie_merge'
include { GFFREAD as GFFREAD_HYBRID_TRANSCRIPTOME           } from '../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HYBRID } from '../../modules/nf-core/star/genomegenerate'
include { GTF_HYBRIDMERGE_GFFCOMPARE } from '../../subworkflows/nf-core/gtf_hybridmerge_gffcompare/main'
include { FASTA_GTF_BAM_RPBP       } from '../../subworkflows/nf-core/fasta_gtf_bam_rpbp/main'
include { GEDI_INDEXGENOME         } from '../../modules/nf-core/gedi/indexgenome/main'
include { GEDI_PRICE               } from '../../modules/nf-core/gedi/price/main'
include { ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE } from '../../subworkflows/nf-core/orftable_fasta_gtf_buildorfcatalogue/main'
include { QUANTIFY_ORF_PSITE       } from '../../subworkflows/local/quantify_orf_psite'
include { DTE_COUNTS_PREP          } from '../../modules/local/dte_counts_prep'
include { DESEQ2_DELTATE as DESEQ2_DELTATE_ORF } from '../../modules/local/deseq2/deltate'
include { ORF_TO_GENE_CDS_COUNTS   } from '../../modules/local/orf_to_gene_cds_counts'
include { GAWK as FILTER_COUNTS_CANONICAL                      } from '../../modules/nf-core/gawk'

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
    // Branch BAMs by sample type so the StringTie path can prefer RNA-seq
    // when available, and so downstream blocks can route Ribo-seq vs RNA-seq
    // independently.
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

    //
    // SUBWORKFLOW: Novel transcript discovery and hybrid GTF construction.
    //
    // Two upstream sources feed the same downstream chain
    // (gffcompare -> class-code filter -> optional rRNA blacklist intersect
    //  -> concat with canonical):
    //   1. `--novel_gtf <file>`  : user-supplied annotation, skip StringTie.
    //   2. `--skip_stringtie false`: run StringTie on RNA-seq BAMs (preferred)
    //      or fall back to Ribo-seq BAMs with tightened parameters.
    //
    // When neither path produces novel transcripts, ch_hybrid_gtf falls back
    // to ch_canonical_gtf so downstream wiring stays uniform.
    //

    ch_hybrid_gtf = ch_canonical_gtf

    def run_stringtie = !params.skip_stringtie && !params.novel_gtf
    def has_user_novel_gtf = params.novel_gtf as Boolean

    if (run_stringtie || has_user_novel_gtf) {

        ch_novel_pre_filter = Channel.empty()

        if (has_user_novel_gtf) {
            ch_novel_pre_filter = Channel
                .fromPath(params.novel_gtf, checkIfExists: true)
                .map { gtf -> [ [id: 'stringtie_merge'], gtf ] }
        }
        else {
            // Prefer RNA-seq BAMs; fall back to Ribo-seq with tightened args.
            // The fallback is signalled via meta.stringtie_fallback so the
            // `withName: STRINGTIE_STRINGTIE` block in conf/modules.config can
            // swap ext.args based on the BAM source.
            def ribo_fallback_args = params.extra_stringtie_args ?: params.stringtie_ribo_fallback_args

            ch_rnaseq_count = ch_genome_bam_by_type.rnaseq.count()

            ch_stringtie_rnaseq = ch_genome_bam_by_type.rnaseq
                .map { meta, bam -> [ meta + [stringtie_fallback: false], bam ] }

            // Ribo-seq is only consumed when no RNA-seq BAMs are present.
            // `combine` with the count value channel gates the emission.
            ch_stringtie_ribo = ch_genome_bam_by_type.riboseq
                .combine(ch_rnaseq_count)
                .filter { _meta, _bam, rna_count -> rna_count == 0 }
                .map { meta, bam, _rna_count ->
                    [ meta + [stringtie_fallback: true], bam ]
                }

            ch_rnaseq_count
                .filter { it == 0 }
                .subscribe {
                    log.warn "No RNA-seq BAMs available for StringTie assembly. Using Ribo-seq BAMs with conservative parameters (${ribo_fallback_args}). Novel transcript assemblies may contain artefacts; review carefully and consider supplying --novel_gtf from a dedicated RNA-seq run."
                }

            BAM_STRINGTIE_MERGE(
                ch_stringtie_rnaseq.mix(ch_stringtie_ribo),
                ch_gtf.map { gtf -> [ [:], gtf ] }
            )

            ch_novel_pre_filter = BAM_STRINGTIE_MERGE.out.stringtie_gtf
                .map { _meta, gtf -> [ [id: 'stringtie_merge'], gtf ] }
        }

        //
        // Classify novel transcripts against the full reference GTF, drop
        // class codes not in --stringtie_class_codes, optionally subtract a
        // strand-aware blacklist, then concatenate with the canonical backbone
        // (synthesising parent gene rows for any novel transcripts whose
        // gene_id isn't already in the backbone). One subworkflow call covers
        // gffcompare + gawk filter + (optional) bedtools intersect + gawk concat.
        //
        ch_blacklist_bed = params.rrna_blacklist
            ? Channel.fromPath(params.rrna_blacklist, checkIfExists: true).map { bed -> [ [id: 'rrna_blacklist'], bed ] }
            : Channel.empty()

        if (!params.rrna_blacklist) {
            log.info "No rRNA/repeat blacklist supplied via --rrna_blacklist; skipping post-assembly blacklist intersect."
        }

        GTF_HYBRIDMERGE_GFFCOMPARE(
            ch_novel_pre_filter,
            ch_gtf.map { gtf -> [ [:], gtf ] }.first(),
            ch_canonical_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] }.first(),
            Channel.value(params.stringtie_class_codes),
            ch_blacklist_bed
        )

        ch_hybrid_gtf = GTF_HYBRIDMERGE_GFFCOMPARE.out.hybrid_gtf.map { _meta, gtf -> gtf }
    }

    //
    // Extended ORF discovery: second STAR pass against a hybrid transcriptome
    // (issue #171). RiboCode requires a transcriptome-coordinate BAM keyed to
    // whichever transcriptome FASTA was used at alignment time. To bring novel
    // intergenic transcripts into RiboCode, we rebuild the transcriptome FASTA
    // from the hybrid GTF and re-align Ribo-seq reads against it.
    //
    // Compute cost: roughly doubles STAR alignment work for Ribo-seq samples.
    // The hybrid transcriptome FASTA and hybrid STAR index are each built once
    // per pipeline run (value channels). The second STAR pass runs only on
    // Ribo-seq samples — RNA-seq and TI-seq are not consumed by RiboCode.
    // riboWaltz stays on the canonical/reference transcriptome BAM by design:
    // it's a QC/calibration tool whose frame plots and metaheatmaps are driven
    // by CDS-bearing canonical transcripts. Feeding it CDS-absent novel
    // transcripts would degrade diagnostic plots without any ORF-discovery gain.
    // See docs/usage.md (Extended ORF discovery) for the full rationale.
    //
    def novel_source_configured = !params.skip_stringtie || params.novel_gtf
    def extended_orf_active = params.extended_orf_analysis && novel_source_configured

    if (params.extended_orf_analysis && !novel_source_configured) {
        log.warn "--extended_orf_analysis is enabled but no novel-transcript source is configured (--skip_stringtie is true and --novel_gtf is unset). The flag has no effect; ORF callers will run against the canonical GTF as usual."
    }

    if (extended_orf_active && params.skip_plastid) {
        log.warn "--extended_orf_analysis is enabled but --skip_plastid is true. ORF-level P-site quantification needs the plastid wiggle tracks and will be skipped; the ORF catalogue will still be built."
    }

    ch_hybrid_transcriptome_bam = Channel.empty()

    if (extended_orf_active) {
        // Extract spliced transcript sequences (canonical + novel) from the
        // hybrid GTF into a single transcriptome FASTA, then rebuild a STAR
        // index against the same genome FASTA + hybrid GTF so a second STAR
        // pass can emit a hybrid transcriptome BAM that RiboCode can consume.
        GFFREAD_HYBRID_TRANSCRIPTOME(
            ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] },
            ch_fasta
        )

        STAR_GENOMEGENERATE_HYBRID(
            ch_fasta.map { fasta -> [ [id: 'hybrid_reference'], fasta ] },
            ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] }
        )

        ch_hybrid_transcriptome_fasta = GFFREAD_HYBRID_TRANSCRIPTOME.out.gffread_fasta.map { _meta, fasta -> fasta }.first()
        ch_hybrid_star_index          = STAR_GENOMEGENERATE_HYBRID.out.index.map { _meta, index -> index }.first()

        // Only Ribo-seq samples feed RiboCode; skip the costly re-alignment
        // for RNA-seq and TI-seq.
        ch_reads_for_hybrid_alignment = ch_reads_for_alignment
            .filter { meta, _reads -> meta.sample_type == 'riboseq' }

        FASTQ_ALIGN_STAR_HYBRID(
            ch_reads_for_hybrid_alignment,
            ch_hybrid_star_index.map { [ [:], it ] },
            ch_hybrid_gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            ch_fasta.map { [ [:], it ] },
            ch_hybrid_transcriptome_fasta.map { [ [:], it ] }
        )

        ch_hybrid_transcriptome_bam = FASTQ_ALIGN_STAR_HYBRID.out.orig_bam_transcript

        ch_multiqc_files = ch_multiqc_files
            .mix(FASTQ_ALIGN_STAR_HYBRID.out.stats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR_HYBRID.out.flagstat.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR_HYBRID.out.idxstats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR_HYBRID.out.log_final.collect{it[1]})
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

    ch_bams_for_analysis = ch_genome_bam_by_type.riboseq.join(ch_genome_bam_index)
    // Use the canonical (one-transcript-per-gene) annotation backbone for ORF
    // calling, P-site calibration and DTE; the full `ch_gtf` is reserved for
    // genome-guided alignment.
    ch_fasta_gtf = ch_fasta.combine(ch_canonical_gtf).map{ fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }.first()
    ch_fasta_gtf_for_ribotish = ch_fasta_gtf.map{ meta, fasta, gtf -> [ meta, fasta, gtf, [] ] }.first()

    //
    // Extended ORF discovery (issues #165 + #171): when --extended_orf_analysis
    // is on and a novel-transcript source is configured, route genome-BAM ORF
    // callers (Ribo-TISH predict, Ribotricer prepare-orfs) to the hybrid GTF
    // so that novel intergenic ORFs are discovered. Ribo-TISH additionally
    // receives the canonical backbone via -a for background + classification.
    // RiboCode additionally consumes the hybrid transcriptome BAM produced by
    // the second STAR pass above (#171); riboWaltz, plastid and Salmon-based
    // quantification stay on the canonical backbone by design.
    // `novel_source_configured` and `extended_orf_active` are defined earlier,
    // alongside the second STAR pass that builds ch_hybrid_transcriptome_bam.
    //
    ch_fasta_gtf_extended = ch_fasta
        .combine(ch_hybrid_gtf)
        .map { fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }
        .first()
    ch_fasta_gtf_for_ribotish_extended = ch_fasta
        .combine(ch_hybrid_gtf)
        .map { fasta, hybrid -> [ [id: 'reference'], fasta, hybrid, [] ] }
        .first()

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

        def ch_ribotish_predict_annotation = extended_orf_active ?
            ch_fasta_gtf_for_ribotish_extended :
            ch_fasta_gtf_for_ribotish

        RIBOTISH_PREDICT_INDIVIDUAL(
            ribotish_predict_inputs.bam,
            [[:],[],[]],
            ch_ribotish_predict_annotation,
            [[:],[]],
            ribotish_predict_inputs.offset,
            [[:],[]]
        )
        ch_versions = ch_versions.mix(RIBOTISH_PREDICT_INDIVIDUAL.out.versions)

        RIBOTISH_PREDICT_ALL(
            ribotish_predict_inputs.bam.map{meta, bam, bai -> [[id:'allsamples'], bam, bai]}.groupTuple(),
            [[:],[],[]],
            ch_ribotish_predict_annotation,
            [[:],[]],
            ribotish_predict_inputs.offset.map{meta, offset -> [[id:'allsamples'], offset]}.groupTuple(),
            [[:],[]]
        )
        ch_versions = ch_versions.mix(RIBOTISH_PREDICT_ALL.out.versions)
    }

    if (params.run_ribotricer){
        log.warn "Ribotricer is enabled via --run_ribotricer. Benchmark data (FK/NGB, May 2026) found its ORF-score column is rank-unstable across biological replicates (mean Spearman 0.288 vs Jaccard 0.770). Its binary calls are usable, but do not rely on its scores as the primary ranking source; the cross-caller rank aggregation will exclude them."

        def ch_ribotricer_annotation = extended_orf_active ?
            ch_fasta_gtf_extended :
            ch_fasta_gtf

        RIBOTRICER_PREPAREORFS(
            ch_ribotricer_annotation
        )
        ch_versions = ch_versions.mix(RIBOTRICER_PREPAREORFS.out.versions)

        RIBOTRICER_DETECTORFS(
            ch_bams_for_analysis,
            RIBOTRICER_PREPAREORFS.out.candidate_orfs
        )
        ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS.out.versions)
    }

    if (params.run_rpbp){
        log.warn "Rp-Bp is enabled via --run_rpbp. Benchmark data (FK/NGB, May 2026) places it at Tier-1 for both rank concordance (mean Spearman 0.893) and set overlap (mean Jaccard 0.673), but expect ~20-24h per replicate at genome-wide scale because the Bayesian MCMC fit dominates. Plan compute accordingly."

        def ch_rpbp_annotation = extended_orf_active ?
            ch_fasta_gtf_extended :
            ch_fasta_gtf

        FASTA_GTF_BAM_RPBP(
            ch_bams_for_analysis,
            ch_rpbp_annotation
        )
    }

    if (params.run_price){
        log.warn "PRICE is enabled via --run_price. PRICE (Erhard et al. 2018) estimates a shared cohort-level codon-position model via EM and is opt-in because its runtime at genome-wide scale is comparable to Rp-Bp. Plan compute accordingly."

        def ch_price_fasta_gtf = extended_orf_active ?
            ch_fasta_gtf_extended :
            ch_fasta_gtf

        GEDI_INDEXGENOME(
            ch_price_fasta_gtf
        )

        // PRICE estimates the codon-position model from the riboseq cohort
        // as a whole, so feed all riboseq BAMs into a single PRICE call.
        ch_price_inputs = ch_bams_for_analysis
            .map { meta, bam, bai -> [ [id: 'allsamples'], bam, bai ] }
            .groupTuple()

        GEDI_PRICE(
            ch_price_inputs,
            GEDI_INDEXGENOME.out.index.first()
        )
    }

    //
    // Dynamic ORF-caller set for cross-caller agreement (issue #07).
    // The enabled list reflects which callers ran at runtime; the agreement
    // threshold and rank-aggregation set are derived from it so the logic
    // works whether 2 (default) or 3+ callers are active.
    //
    def enabled_orf_callers = []
    if (!params.skip_ribotish)   { enabled_orf_callers << 'ribotish' }
    if (!params.skip_ribocode)   { enabled_orf_callers << 'ribocode' }
    if ( params.run_ribotricer)  { enabled_orf_callers << 'ribotricer' }
    if ( params.run_rpbp)        { enabled_orf_callers << 'rpbp' }
    if ( params.run_price)       { enabled_orf_callers << 'price' }

    // Ribotricer contributes binary calls only; its scores are excluded from
    // the cross-caller rank aggregation due to known rank instability. Rp-Bp's
    // Bayes factor is stable (Spearman 0.893) and is retained for ranking.
    def rank_aggregation_callers  = enabled_orf_callers - 'ribotricer'
    // Strict-majority of enabled callers (floor(N/2)+1): N=2 -> 2 (both must
    // agree), N=3 -> 2 (majority). Adapts as the caller set grows.
    def orf_agreement_min_callers = enabled_orf_callers
        ? enabled_orf_callers.size().intdiv(2) + 1
        : 0
    ch_enabled_orf_callers      = Channel.value(enabled_orf_callers)
    ch_rank_aggregation_callers = Channel.value(rank_aggregation_callers)

    if (!params.skip_ribocode){
        // RiboCode requires transcriptome-coordinate BAMs. When extended-ORF
        // analysis is active (#171), swap in the hybrid transcriptome BAM
        // (second STAR pass against the hybrid GTF, Ribo-seq only) plus the
        // hybrid GTF as the annotation source so novel intergenic transcripts
        // are visible to RiboCode. Otherwise keep the canonical wiring.
        ch_ribocode_transcriptome_bam_source = extended_orf_active ?
            ch_hybrid_transcriptome_bam :
            ch_transcriptome_bam

        ch_transcriptome_bams_for_ribocode = ch_ribocode_transcriptome_bam_source
            .branch { meta, bam ->
                riboseq: meta.sample_type == 'riboseq'
                    return [ meta, bam ]
            }
            .riboseq

        def ch_ribocode_gtf_source = extended_orf_active ?
            ch_hybrid_gtf :
            ch_gtf

        // Step 1: Update GTF annotation
        def ribocode_gtf_meta_id = extended_orf_active ? 'hybrid_reference' : 'reference'
        RIBOCODE_GTFUPDATE(
            ch_ribocode_gtf_source.map { [ [id: ribocode_gtf_meta_id], it ] }.first()
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
        ch_ribocode_with_config = ch_transcriptome_bams_for_ribocode
            .join(RIBOCODE_METAPLOTS.out.config, remainder: true)

        ch_ribocode_with_config
            .filter { _meta, _bam, config -> config == null }
            .subscribe { meta, _bam, _config ->
                log.warn "Skipping ${meta.id} for RiboCode: RIBOCODE_METAPLOTS produced no config (likely sparse periodicity data)."
            }

        ch_ribocode_inputs = ch_ribocode_with_config
            .filter { _meta, _bam, config -> config != null }

        RIBOCODE_RIBOCODE(
            ch_ribocode_inputs.map { meta, bam, _config  -> [ meta, bam ] },
            RIBOCODE_PREPARE.out.annotation,
            ch_ribocode_inputs.map { meta, _bam, config  -> [ meta, config ] }
        )
    }

    //
    // Cross-sample ORF catalogue (issue #167). Built once per pipeline run
    // when extended-ORF analysis is enabled AND at least one ORF caller ran.
    // The catalogue normalises each caller's per-sample output into a unified
    // BED12, merges with class-aware collapse (transcript-ID for canonical
    // CDS, reciprocal overlap for novel intergenic / smORFs), and emits the
    // merged BED12 + sidecar TSV + ORF-to-gene mapping + AA FASTA.
    //
    if (extended_orf_active && enabled_orf_callers) {
        ch_ribotish_pred_cat   = (!params.skip_ribotish)  ? RIBOTISH_PREDICT_INDIVIDUAL.out.predictions : Channel.empty()
        ch_ribocode_pred_cat   = (!params.skip_ribocode)  ? RIBOCODE_RIBOCODE.out.orf_txt              : Channel.empty()
        ch_ribotricer_pred_cat = ( params.run_ribotricer) ? RIBOTRICER_DETECTORFS.out.orfs             : Channel.empty()
        ch_rpbp_pred_cat       = ( params.run_rpbp)       ? FASTA_GTF_BAM_RPBP.out.predicted           : Channel.empty()
        ch_price_pred_cat      = ( params.run_price)      ? GEDI_PRICE.out.orfs_tsv                   : Channel.empty()

        ch_orf_tables = ch_ribotish_pred_cat  .map { meta, f -> [ meta, f, 'ribotish'   ] }
            .mix(ch_ribocode_pred_cat         .map { meta, f -> [ meta, f, 'ribocode'   ] })
            .mix(ch_ribotricer_pred_cat       .map { meta, f -> [ meta, f, 'ribotricer' ] })
            .mix(ch_rpbp_pred_cat             .map { meta, f -> [ meta, f, 'rpbp'       ] })
            .mix(ch_price_pred_cat            .map { meta, f -> [ meta, f, 'price'      ] })

        ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE(
            ch_orf_tables,
            ch_fasta    .map { fasta -> [ [id: 'reference'], fasta ] }.first(),
            ch_hybrid_gtf.map { gtf   -> [ [id: 'reference'], gtf   ] }.first()
        )
        ch_multiqc_files = ch_multiqc_files.mix(ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.multiqc.collect{it[1]}.ifEmpty([]))
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

    // ORF-level count matrix; populated below when extended ORF analysis is
    // active and plastid is enabled. Consumed by issue #168 (DTE).
    ch_orf_count_matrix = Channel.empty()

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

        //
        // Per-ORF P-site quantification (issue #166). Runs additively to the
        // gene-level QUANTIFY_INFRAME_PSITE_PLASTID path. Gated on the same
        // predicate as ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE so it only fires when the
        // catalogue exists.
        //
        if (extended_orf_active && enabled_orf_callers) {
            ch_orf_psite_tracks = PLASTID_MAKE_WIGGLE.out.tracks
                .map { meta, tracks -> [ meta, tracks[0], tracks[1] ] }

            QUANTIFY_ORF_PSITE (
                ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.catalogue_bed12,
                ch_orf_psite_tracks,
                ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.orf_to_gene_tsv
            )
            ch_versions         = ch_versions.mix(QUANTIFY_ORF_PSITE.out.versions)
            ch_orf_count_matrix = QUANTIFY_ORF_PSITE.out.orf_count_matrix
        }

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
            file("${projectDir}/assets/gtf_to_inframe_psites.awk"),
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

        // Issue #168 Tier 1: when extended ORF analysis is active and a
        // cohort ORF P-site matrix exists, replace the plastid-derived
        // gene-CDS p-site counts with a re-aggregation that sums ONLY
        // canonical_cds ORFs from the catalogue. This keeps the gene-level
        // TE numerator clean of uORF/dORF dynamics, which are picked up
        // separately in the Tier 2 ORF-level DTE below.
        if (extended_orf_active && enabled_orf_callers) {
            ORF_TO_GENE_CDS_COUNTS(
                ch_orf_count_matrix
                    .combine(ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.orf_to_gene_tsv.map { _meta, tsv -> tsv })
                    .combine(ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.catalogue_tsv.map { _meta, tsv -> tsv })
                    .map { meta, orf_counts, o2g, cat_tsv -> [meta, orf_counts, o2g, cat_tsv] }
            )
            ch_psite_counts_merged = ORF_TO_GENE_CDS_COUNTS.out.gene_counts
        }

        REPLACE_RIBOSEQ_COUNTS_IN_MATRIX(
            ch_psite_counts_merged
                .combine(QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map{ meta, counts -> counts })
                .map { meta, psite_counts, salmon_counts -> [meta, [psite_counts, salmon_counts]] },
            file("${projectDir}/assets/replace_riboseq_counts_in_matrix.awk"),
            false
        )
        ch_te_counts = REPLACE_RIBOSEQ_COUNTS_IN_MATRIX.out.output
        ch_versions = ch_versions.mix(QUANTIFY_INFRAME_PSITE_PLASTID.out.versions)
    }

    //
    // Do a translational efficiency analysis where contrasts are supplied
    //

    if (ch_contrasts_file){

        // GTF first, counts second: the awk script uses the FNR==NR idiom
        // to build the gene_id keep-set from the first file (the GTF) before
        // streaming the second (the counts TSV) and emitting matched rows.
        ch_filter_counts_in = ch_te_counts
            .combine(ch_canonical_gtf)
            .map { meta, counts, gtf -> [ meta, [ gtf, counts ] ] }

        FILTER_COUNTS_CANONICAL(
            ch_filter_counts_in,
            file("${projectDir}/assets/filter_counts_canonical.awk"),
            false
        )
        ch_te_counts = FILTER_COUNTS_CANONICAL.out.output

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

        if (extended_orf_active && enabled_orf_callers && !params.skip_plastid) {
            DTE_COUNTS_PREP(
                ch_orf_count_matrix
                    .combine(QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map { _meta, counts -> counts })
                    .combine(ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE.out.orf_to_gene_tsv.map { _meta, tsv -> tsv })
                    .map { meta, ribo, rna, o2g -> [ meta, ribo, rna, o2g ] }
            )

            ch_orf_samplesheet_matrix = DTE_COUNTS_PREP.out.counts
                .combine(ch_samplesheet)
                .map { meta, counts, samplesheet -> [ meta, samplesheet, counts ] }
                .first()

            DESEQ2_DELTATE_ORF(
                ch_contrasts,
                ch_orf_samplesheet_matrix
            )
            ch_versions = ch_versions.mix(DTE_COUNTS_PREP.out.versions)
            ch_versions = ch_versions.mix(DESEQ2_DELTATE_ORF.out.versions)
        }
    }

    // Issue #168 Tier 3 (DOTSeq): deferred. The package is in
    // Bioconductor devel only; the param is a placeholder so users can
    // express intent now and so the schema/docs stay aligned with the
    // implementation plan. Tracked in #168.
    if (params.run_dotseq) {
        log.info "--run_dotseq is enabled, but DOTSeq Tier-3 ORF-level DTE/DOU is deferred to a future release (DOTSeq is currently in Bioconductor devel only). See issue #168."
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
    multiqc_report   = ch_multiqc_report   // channel: /path/to/multiqc_report.html
    versions         = ch_versions         // channel: [ path(versions.yml) ]
    hybrid_gtf       = ch_hybrid_gtf       // channel: path(hybrid_reference.gtf) - canonical + filtered novel; equals canonical when no novel source is configured
    orf_count_matrix = ch_orf_count_matrix // channel: [ meta, orf_psite_counts.tsv ] - per-ORF P-site count matrix (issue #166); empty unless extended-ORF + plastid both active
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
