//
// Build a cohort-level ORF catalogue from per-sample, per-caller outputs.
//
// Chains:
//   1. Per-caller normalisation (BED12 + sidecar TSV).
//   2. Concatenation across samples + callers into a single (BED, TSV) pair.
//   3. ORF_CATALOGUE_MERGE applies class-aware collapse strategy.
//   4. BEDTOOLS_GETFASTA + SEQKIT_TRANSLATE produce the AA FASTA from the
//      merged BED12 (-split -s -nameOnly walks BED12 blocks in mRNA order
//      and reverse-complements minus-strand records; seqkit translate --trim
//      drops trailing stop codons).
//
// Inputs that come from optional callers (ribotricer, rpbp, price) may be
// passed as `Channel.empty()`; the subworkflow handles missing callers cleanly.
//
include { ORF_NORMALISE as ORF_NORMALISE_RIBOTISH    } from '../../../modules/local/orf_normalise'
include { ORF_NORMALISE as ORF_NORMALISE_RIBOCODE    } from '../../../modules/local/orf_normalise'
include { ORF_NORMALISE as ORF_NORMALISE_RIBOTRICER  } from '../../../modules/local/orf_normalise'
include { ORF_NORMALISE as ORF_NORMALISE_RPBP        } from '../../../modules/local/orf_normalise'
include { ORF_NORMALISE as ORF_NORMALISE_PRICE       } from '../../../modules/local/orf_normalise'
include { ORF_CATALOGUE_MERGE                        } from '../../../modules/local/orf_catalogue_merge'
include { BEDTOOLS_GETFASTA                          } from '../../../modules/nf-core/bedtools/getfasta/main'
include { SEQKIT_TRANSLATE                           } from '../../../modules/nf-core/seqkit/translate/main'

workflow BUILD_ORF_CATALOGUE {

    take:
    ch_ribotish_pred    // channel: [ val(meta), path(*_pred.txt) ]    - per sample
    ch_ribocode_pred    // channel: [ val(meta), path(*.txt) ]         - per sample
    ch_ribotricer_pred  // channel: [ val(meta), path(*_translating_ORFs.tsv) ] - per sample, may be empty
    ch_rpbp_pred        // channel: [ val(meta), path(*.predicted-orfs.bed.gz) ] - per sample, may be empty
    ch_price_pred       // channel: [ val(meta), path(*.orfs.tsv) ]    - cohort-level (PRICE is one-shot), may be empty
    ch_fasta            // channel: path(genome.fasta)
    ch_gtf              // channel: path(annotation.gtf) - canonical or hybrid

    main:

    // Reference channel for normalisers: [ [id: 'reference'], gtf ]
    ch_ref_gtf = ch_gtf.map { gtf -> [ [id: 'reference'], gtf ] }.first()

    // Tag each input with meta.caller so the consolidated ORF_NORMALISE
    // template can dispatch on caller-specific parsing / classification.
    ORF_NORMALISE_RIBOTISH   ( ch_ribotish_pred.map   { meta, f -> [ meta + [caller: 'ribotish'],   f ] }, ch_ref_gtf )
    ORF_NORMALISE_RIBOCODE   ( ch_ribocode_pred.map   { meta, f -> [ meta + [caller: 'ribocode'],   f ] }, ch_ref_gtf )
    ORF_NORMALISE_RIBOTRICER ( ch_ribotricer_pred.map { meta, f -> [ meta + [caller: 'ribotricer'], f ] }, ch_ref_gtf )
    ORF_NORMALISE_RPBP       ( ch_rpbp_pred.map       { meta, f -> [ meta + [caller: 'rpbp'],       f ] }, ch_ref_gtf )
    ORF_NORMALISE_PRICE      ( ch_price_pred.map      { meta, f -> [ meta + [caller: 'price'],      f ] }, ch_ref_gtf )

    // Collect normalised BED12s and sidecar TSVs across callers and samples
    // into a single cohort-keyed channel: [ [id:'allsamples'], [beds...], [tsvs...] ].
    ch_all_beds = ORF_NORMALISE_RIBOTISH.out.bed12
        .mix(ORF_NORMALISE_RIBOCODE.out.bed12)
        .mix(ORF_NORMALISE_RIBOTRICER.out.bed12)
        .mix(ORF_NORMALISE_RPBP.out.bed12)
        .mix(ORF_NORMALISE_PRICE.out.bed12)
        .map { _meta, bed -> bed }
        .collect()
        .map { beds -> [ 'allsamples', beds ] }

    ch_all_tsvs = ORF_NORMALISE_RIBOTISH.out.tsv
        .mix(ORF_NORMALISE_RIBOCODE.out.tsv)
        .mix(ORF_NORMALISE_RIBOTRICER.out.tsv)
        .mix(ORF_NORMALISE_RPBP.out.tsv)
        .mix(ORF_NORMALISE_PRICE.out.tsv)
        .map { _meta, tsv -> tsv }
        .collect()
        .map { tsvs -> [ 'allsamples', tsvs ] }

    ch_merge_in = ch_all_beds
        .combine(ch_all_tsvs, by: 0)
        .map { _key, beds, tsvs -> [ [id: 'allsamples'], beds, tsvs ] }

    ORF_CATALOGUE_MERGE ( ch_merge_in )

    BEDTOOLS_GETFASTA (
        ORF_CATALOGUE_MERGE.out.catalogue_bed,
        ch_fasta.first()
    )

    SEQKIT_TRANSLATE ( BEDTOOLS_GETFASTA.out.fasta )

    emit:
    catalogue_bed    = ORF_CATALOGUE_MERGE.out.catalogue_bed
    catalogue_tsv    = ORF_CATALOGUE_MERGE.out.catalogue_tsv
    catalogue_faa    = SEQKIT_TRANSLATE.out.fastx
    orf_to_gene_tsv  = ORF_CATALOGUE_MERGE.out.orf_to_gene
    multiqc_summary  = ORF_CATALOGUE_MERGE.out.mqc
}
