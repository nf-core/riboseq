//
// Build a cohort-level ORF catalogue from per-sample, per-caller outputs.
//
// Chains:
//   1. Per-caller normalisation (BED12 + sidecar TSV).
//   2. Concatenation across samples + callers into a single (BED, TSV) pair.
//   3. ORF_CATALOGUE_MERGE applies class-aware collapse strategy.
//   4. ORF_CATALOGUE_EXTRACTAA produces the AA FASTA from the merged BED12.
//
// Inputs that come from optional callers (ribotricer, rpbp) may be passed
// as `Channel.empty()`; the subworkflow handles missing callers cleanly.
//
include { ORF_NORMALISE_RIBOTISH    } from '../../../modules/local/orf/normalise/ribotish'
include { ORF_NORMALISE_RIBOCODE    } from '../../../modules/local/orf/normalise/ribocode'
include { ORF_NORMALISE_RIBOTRICER  } from '../../../modules/local/orf/normalise/ribotricer'
include { ORF_NORMALISE_RPBP        } from '../../../modules/local/orf/normalise/rpbp'
include { ORF_CATALOGUE_MERGE       } from '../../../modules/local/orf/catalogue/merge'
include { ORF_CATALOGUE_EXTRACTAA  } from '../../../modules/local/orf/catalogue/extractaa'

workflow BUILD_ORF_CATALOGUE {

    take:
    ch_ribotish_pred    // channel: [ val(meta), path(*_pred.txt) ]    - per sample
    ch_ribocode_pred    // channel: [ val(meta), path(*.txt) ]         - per sample
    ch_ribotricer_pred  // channel: [ val(meta), path(*_translating_ORFs.tsv) ] - per sample, may be empty
    ch_rpbp_pred        // channel: [ val(meta), path(*.predicted-orfs.bed.gz) ] - per sample, may be empty
    ch_fasta            // channel: path(genome.fasta)
    ch_gtf              // channel: path(annotation.gtf) - canonical or hybrid

    main:

    ch_versions = Channel.empty()

    // Reference channel for normalisers: [ [id: 'reference'], gtf ]
    ch_ref_gtf = ch_gtf.map { gtf -> [ [id: 'reference'], gtf ] }.first()

    ORF_NORMALISE_RIBOTISH ( ch_ribotish_pred,   ch_ref_gtf )
    ORF_NORMALISE_RIBOCODE ( ch_ribocode_pred,   ch_ref_gtf )
    ORF_NORMALISE_RIBOTRICER ( ch_ribotricer_pred, ch_ref_gtf )
    ORF_NORMALISE_RPBP     ( ch_rpbp_pred,       ch_ref_gtf )

    // Collect normalised BED12s and sidecar TSVs across callers and samples
    // into a single cohort-keyed channel: [ [id:'allsamples'], [beds...], [tsvs...] ].
    ch_all_beds = ORF_NORMALISE_RIBOTISH.out.bed12
        .mix(ORF_NORMALISE_RIBOCODE.out.bed12)
        .mix(ORF_NORMALISE_RIBOTRICER.out.bed12)
        .mix(ORF_NORMALISE_RPBP.out.bed12)
        .map { _meta, bed -> bed }
        .collect()
        .ifEmpty([])

    ch_all_tsvs = ORF_NORMALISE_RIBOTISH.out.tsv
        .mix(ORF_NORMALISE_RIBOCODE.out.tsv)
        .mix(ORF_NORMALISE_RIBOTRICER.out.tsv)
        .mix(ORF_NORMALISE_RPBP.out.tsv)
        .map { _meta, tsv -> tsv }
        .collect()
        .ifEmpty([])

    ch_merge_in = ch_all_beds
        .combine(ch_all_tsvs)
        .map { beds, tsvs -> [ [id: 'allsamples'], beds, tsvs ] }

    ORF_CATALOGUE_MERGE ( ch_merge_in )

    ORF_CATALOGUE_EXTRACTAA (
        ORF_CATALOGUE_MERGE.out.catalogue_bed,
        ch_fasta.map { fasta -> [ [id: 'reference'], fasta ] }.first()
    )

    emit:
    catalogue_bed    = ORF_CATALOGUE_MERGE.out.catalogue_bed
    catalogue_tsv    = ORF_CATALOGUE_MERGE.out.catalogue_tsv
    catalogue_faa    = ORF_CATALOGUE_EXTRACTAA.out.faa
    orf_to_gene_tsv  = ORF_CATALOGUE_MERGE.out.orf_to_gene
    multiqc_summary  = ORF_CATALOGUE_MERGE.out.mqc
    versions         = ch_versions
}
