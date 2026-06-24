//
// Conditional ORF-caller dispatch for the riboseq pipeline.
//
// Runs each enabled caller against the appropriate annotation: when extended
// ORF analysis is active (issue #165), the genome-BAM callers (Ribo-TISH
// predict, Ribotricer) receive the hybrid GTF so that novel intergenic ORFs
// are within scope. Otherwise everything stays on the canonical backbone.
// Ribo-TISH additionally takes the canonical backbone via -a for background +
// classification. RiboCode is a transcriptome-coordinate caller and stays on
// the canonical wiring here.
//
// Per-caller gating (params.skip_ribotish / params.run_ribotricer) lives here.
//
// Emits one prediction channel per caller plus collected versions and
// multiqc files. Predictions are empty channels for callers that did not run.
//

include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_RIBOSEQ    } from '../../../modules/nf-core/ribotish/quality'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_INDIVIDUAL } from '../../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_ALL        } from '../../../modules/nf-core/ribotish/predict'
include { RIBOTRICER_PREPAREORFS                          } from '../../../modules/nf-core/ribotricer/prepareorfs'
include { RIBOTRICER_DETECTORFS                           } from '../../../modules/nf-core/ribotricer/detectorfs'
include { RIBOCODE_GTFUPDATE                              } from '../../../modules/nf-core/ribocode/gtfupdate'
include { RIBOCODE_PREPARE                                } from '../../../modules/nf-core/ribocode/prepare'
include { RIBOCODE_METAPLOTS                              } from '../../../modules/nf-core/ribocode/metaplots'
include { RIBOCODE_RIBOCODE                               } from '../../../modules/nf-core/ribocode/ribocode'

workflow ORF_CALLER_DISPATCH {

    take:
    ch_bams_for_analysis     // channel: [ val(meta), path(bam), path(bai) ] - riboseq genome BAMs
    ch_transcriptome_bam     // channel: [ val(meta), path(bam) ] - all sample types, canonical transcriptome
    ch_fasta                 // channel: path(genome.fasta)
    ch_canonical_gtf         // channel: path(canonical.gtf)
    ch_hybrid_gtf            // channel: path(hybrid.gtf)        - equals canonical when no novel source
    ch_gtf                   // channel: path(genome.gtf)        - full multi-isoform, used by RiboCode
    extended_orf_active      // val: bool

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Annotation channels. Canonical for ORF calling / P-site / DTE; the full
    // ch_gtf is reserved for genome-guided alignment elsewhere in the pipeline.
    ch_fasta_gtf = ch_fasta.combine(ch_canonical_gtf).map{ fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }.first()
    ch_fasta_gtf_for_ribotish = ch_fasta_gtf.map{ meta, fasta, gtf -> [ meta, fasta, gtf, [] ] }.first()

    ch_fasta_gtf_extended = ch_fasta
        .combine(ch_hybrid_gtf)
        .map { fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }
        .first()
    // Extended mode: hybrid GTF feeds Ribo-TISH `-g` (discovery target), canonical
    // backbone feeds `-a` (background model + ORF classification labels).
    ch_fasta_gtf_for_ribotish_extended = ch_fasta
        .combine(ch_hybrid_gtf)
        .combine(ch_canonical_gtf)
        .map { fasta, hybrid, canonical -> [ [id: 'reference'], fasta, hybrid, canonical ] }
        .first()

    //
    // Ribo-TISH
    //
    ch_ribotish_predictions = Channel.empty()
    if (!params.skip_ribotish) {
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

        // `-a` secondary annotation: canonical backbone in extended mode, empty
        // ([]) otherwise. RIBOTISH_PREDICT guards on a truthy reference_gtf
        // (reference_gtf ? "-a ..." : ''), so an empty payload runs without -a.
        def ch_predict_fasta_gtf = ch_ribotish_predict_annotation.map { meta, fasta, gtf, _ref -> [ meta, fasta, gtf ] }
        def ch_predict_ref_gtf   = ch_ribotish_predict_annotation.map { _meta, _fasta, _gtf, ref -> [ [id: 'reference'], ref ] }

        RIBOTISH_PREDICT_INDIVIDUAL(
            ribotish_predict_inputs.bam,
            [[:],[],[]],
            ch_predict_fasta_gtf,
            [[:],[]],
            ribotish_predict_inputs.offset,
            [[:],[]],
            ch_predict_ref_gtf
        )
        ch_ribotish_predictions = RIBOTISH_PREDICT_INDIVIDUAL.out.predictions

        RIBOTISH_PREDICT_ALL(
            ribotish_predict_inputs.bam.map{ _meta, bam, bai -> [[id:'allsamples'], bam, bai]}.groupTuple(),
            [[:],[],[]],
            ch_predict_fasta_gtf,
            [[:],[]],
            ribotish_predict_inputs.offset.map{ _meta, offset -> [[id:'allsamples'], offset]}.groupTuple(),
            [[:],[]],
            ch_predict_ref_gtf
        )
    }

    //
    // Ribotricer
    //
    ch_ribotricer_predictions = Channel.empty()
    if (params.run_ribotricer) {
        log.warn "Ribotricer is enabled via --run_ribotricer. Its per-ORF scores are unstable across biological replicates, so its binary calls contribute to cross-caller agreement but its scores are excluded from the rank aggregation."

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
        ch_ribotricer_predictions = RIBOTRICER_DETECTORFS.out.orfs
    }

    //
    // RiboCode
    //
    ch_ribocode_predictions = Channel.empty()
    if (!params.skip_ribocode) {
        // RiboCode requires transcriptome-coordinate BAMs keyed to the
        // reference transcriptome FASTA used at alignment time, so it stays on
        // the reference-transcriptome wiring regardless of --extended_orf_analysis.
        ch_transcriptome_bams_for_ribocode = ch_transcriptome_bam
            .branch { meta, bam ->
                riboseq: meta.sample_type == 'riboseq'
                    return [ meta, bam ]
            }
            .riboseq

        // Step 1: Update GTF annotation
        RIBOCODE_GTFUPDATE(
            ch_gtf.map { [ [id: 'reference'], it ] }.first()
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
        ch_ribocode_inputs = ch_transcriptome_bams_for_ribocode
            .join(RIBOCODE_METAPLOTS.out.config)

        RIBOCODE_RIBOCODE(
            ch_ribocode_inputs.map { meta, bam, _config -> [ meta, bam ] },
            RIBOCODE_PREPARE.out.annotation,
            ch_ribocode_inputs.map { meta, _bam, config -> [ meta, config ] }
        )
        ch_ribocode_predictions = RIBOCODE_RIBOCODE.out.orf_txt
    }

    emit:
    ribotish_predictions   = ch_ribotish_predictions   // [ meta, predictions.txt ] or empty
    ribocode_predictions   = ch_ribocode_predictions   // [ meta, orf.txt          ] or empty
    ribotricer_predictions = ch_ribotricer_predictions // [ meta, orfs             ] or empty
    multiqc_files          = ch_multiqc_files
    versions               = ch_versions
}
