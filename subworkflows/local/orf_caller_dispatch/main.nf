//
// Conditional ORF-caller dispatch for the riboseq pipeline.
//
// Runs each enabled caller against the appropriate annotation: when extended
// ORF analysis is active, genome-BAM callers (Ribo-TISH
// predict, Ribotricer, Rp-Bp, PRICE) receive the hybrid GTF and RiboCode
// receives the hybrid transcriptome BAM + hybrid GTF. Otherwise everything
// stays on the canonical backbone. Ribo-TISH additionally takes the canonical
// backbone via -a for background + classification.
//
// Per-caller gating (params.skip_ribotish / params.run_*) lives here; the
// downstream catalogue gating still lives at the call site.
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
include { FASTA_GTF_BAM_RPBP                              } from '../../nf-core/fasta_gtf_bam_rpbp/main'
include { GEDI_INDEXGENOME                                } from '../../../modules/nf-core/gedi/indexgenome/main'
include { GEDI_PRICE                                      } from '../../../modules/nf-core/gedi/price/main'

workflow ORF_CALLER_DISPATCH {

    take:
    ch_bams_for_analysis     // channel: [ val(meta), path(bam), path(bai) ] - riboseq genome BAMs
    ch_transcriptome_bam     // channel: [ val(meta), path(bam) ] - all sample types, canonical transcriptome
    ch_hybrid_transcriptome_bam // channel: [ val(meta), path(bam) ] - riboseq only, hybrid transcriptome (or empty)
    // The four reference channels below must be value channels: they are paired with the per-sample
    // BAM channels in the caller fan-outs, and a single-item queue would serve only one sample.
    ch_fasta                 // value channel: path(genome.fasta)
    ch_canonical_gtf         // value channel: path(canonical.gtf)
    ch_hybrid_gtf            // value channel: path(hybrid.gtf)  - equals canonical when no novel source
    ch_gtf                   // value channel: path(genome.gtf)  - full multi-isoform, used by RiboCode when not extended
    // Same contract; populated whenever the hybrid GTF is built.
    ch_full_hybrid_gtf       // value channel: path(full_hybrid.gtf) - ch_gtf plus the novel genes
    extended_orf_active      // val: bool

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Annotation channels. Canonical for ORF calling / P-site / DTE; the full
    // ch_gtf is reserved for genome-guided alignment elsewhere in the pipeline.
    ch_fasta_gtf = ch_fasta.combine(ch_canonical_gtf).map{ fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }.first()
    ch_fasta_gtf_for_ribotish = ch_fasta_gtf.map{ meta, fasta, gtf -> [ meta, fasta, gtf, [] ] }

    ch_fasta_gtf_extended = ch_fasta
        .combine(ch_hybrid_gtf)
        .map { fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }
        .first()

    // Annotation for the callers that resolve overlapping ORFs themselves.
    ch_fasta_gtf_multi_isoform = extended_orf_active ?
        ch_fasta.combine(ch_full_hybrid_gtf).map { fasta, gtf -> [ [id: 'full_hybrid_reference'], fasta, gtf ] }.first() :
        ch_fasta.combine(ch_gtf).map { fasta, gtf -> [ [id: 'reference'], fasta, gtf ] }.first()

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
            ch_canonical_gtf.map { [ [:], it ] }
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
    // Rp-Bp
    //
    ch_rpbp_predictions = Channel.empty()
    if (params.run_rpbp) {
        log.warn "Rp-Bp is enabled via --run_rpbp. Expect roughly 20-24h per replicate at genome-wide scale because the Bayesian MCMC fit dominates; plan compute accordingly. Its score column (Bayes factor) is stable and is retained in the cross-caller rank aggregation."

        // Rp-Bp enumerates candidate ORFs per transcript isoform and resolves
        // redundancy/overlaps itself (longest-per-stop, best Bayes factor per
        // overlap; Malone et al. 2017), so it takes the multi-isoform annotation
        // rather than the one-transcript-per-gene backbone, which exists to
        // disambiguate P-site quantification. Restricting it to canonical would
        // drop isoform-specific ORFs and bias ORF-type classification to
        // canonical CDS.
        FASTA_GTF_BAM_RPBP(
            ch_bams_for_analysis,
            ch_fasta_gtf_multi_isoform
        )
        ch_rpbp_predictions = FASTA_GTF_BAM_RPBP.out.predicted
    }

    //
    // PRICE
    //
    ch_price_predictions = Channel.empty()
    if (params.run_price) {
        log.warn "PRICE is enabled via --run_price. PRICE (Erhard et al. 2018) estimates a shared cohort-level codon-position model via EM and is opt-in because its genome-wide runtime is substantial. Plan compute accordingly."

        // PRICE resolves overlapping ORFs and rescues multimappers with its own
        // EM (Erhard et al. 2018), so it has no need for the one-transcript-per-gene
        // backbone that exists to disambiguate P-site quantification. Restricting it
        // would only narrow ORF discovery and bias ORF-type classification to
        // canonical CDS.
        GEDI_INDEXGENOME(
            ch_fasta_gtf_multi_isoform
        )

        // PRICE estimates the codon-position model from the riboseq cohort
        // as a whole, so feed all riboseq BAMs into a single PRICE call.
        ch_price_inputs = ch_bams_for_analysis
            .map { _meta, bam, bai -> [ [id: 'allsamples'], bam, bai ] }
            .groupTuple()

        GEDI_PRICE(
            ch_price_inputs,
            GEDI_INDEXGENOME.out.index
        )
        ch_price_predictions = GEDI_PRICE.out.orfs_tsv
    }

    //
    // RiboCode
    //
    ch_ribocode_predictions = Channel.empty()
    if (!params.skip_ribocode) {
        // RiboCode requires transcriptome-coordinate BAMs. When extended-ORF
        // analysis is active, swap in the hybrid transcriptome BAM
        // (second STAR pass against the hybrid GTF, Ribo-seq only) plus the
        // hybrid GTF as the annotation source so novel intergenic transcripts
        // are visible to RiboCode. Otherwise keep the reference-transcriptome wiring.
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
            ch_ribocode_gtf_source.map { [ [id: ribocode_gtf_meta_id], it ] }
        )

        // Step 2: Prepare annotation files
        RIBOCODE_PREPARE(
            ch_fasta.map { [ [:], it ] },
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
    rpbp_predictions       = ch_rpbp_predictions       // [ meta, predicted        ] or empty
    price_predictions      = ch_price_predictions      // [ meta, orfs.tsv         ] or empty
    multiqc_files          = ch_multiqc_files
    versions               = ch_versions
}
