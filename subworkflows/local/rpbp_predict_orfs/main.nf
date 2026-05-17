//
// Rp-Bp post-alignment ORF prediction chain (one Nextflow process per rpbp tool).
//
// We replicate the chain that rpbp's `create-orf-profiles` + `predict-translated-orfs`
// wrappers would run internally - but as separate Nextflow processes so each step
// caches independently and we can re-run only what changed.
//
// rpbp's own wrappers (`run-rpbp-pipeline`, `create-orf-profiles`,
// `predict-translated-orfs`) chain these tools internally and additionally invoke
// flexbar + bowtie + STAR for raw FASTQ input. We do alignment with the standard
// nf-core STAR module and skip those steps entirely.
//
// Order (after RPBP_BUILDCONFIG + RPBP_PREPAREGENOME, both once per run):
//   1. extract-metagene-profiles            -> *.metagene-profile.csv.gz
//   2. estimate-metagene-profile-bayes-factors -> *.metagene-periodicity-bayes-factors.csv.gz
//   3. select-periodic-offsets              -> *.periodic-offsets.csv.gz
//   4. extract-orf-profiles                 -> *.profiles.mtx.gz
//      (filter logic replicated from rpbp's ribo_utils.utils.get_periodic_lengths_and_offsets
//       lives in the EXTRACTORFPROFILES module's script block)
//   5. estimate-orf-bayes-factors           -> *.bayes-factors.bed.gz
//   6. select-final-prediction-set          -> *.predicted-orfs.filtered.{bed.gz,dna.fa,protein.fa}
//

include { RPBP_EXTRACTMETAGENEPROFILES        } from '../../../modules/local/rpbp/extractmetageneprofiles'
include { RPBP_ESTIMATEMETAGENEBAYESFACTORS   } from '../../../modules/local/rpbp/estimatemetagenebayesfactors'
include { RPBP_SELECTPERIODICOFFSETS          } from '../../../modules/local/rpbp/selectperiodicoffsets'
include { RPBP_EXTRACTORFPROFILES             } from '../../../modules/local/rpbp/extractorfprofiles'
include { RPBP_ESTIMATEORFBAYESFACTORS        } from '../../../modules/local/rpbp/estimateorfbayesfactors'
include { RPBP_SELECTFINALPREDICTIONSET       } from '../../../modules/local/rpbp/selectfinalpredictionset'

workflow RPBP_PREDICT_ORFS {

    take:
    ch_bam               // channel: [ val(meta), path(bam), path(bai) ] - per Ribo-seq sample
    ch_transcript_bed    // channel: path - annotated transcripts BED from PREPAREGENOME
    ch_orfs_genomic_bed  // channel: path - ORF genomic BED from PREPAREGENOME
    ch_orfs_exons_bed    // channel: path - ORF exons BED from PREPAREGENOME
    ch_genome_fasta      // channel: path - genome FASTA

    main:

    ch_versions = Channel.empty()

    RPBP_EXTRACTMETAGENEPROFILES (
        ch_bam,
        ch_transcript_bed
    )

    RPBP_ESTIMATEMETAGENEBAYESFACTORS (
        RPBP_EXTRACTMETAGENEPROFILES.out.metagene
    )

    RPBP_SELECTPERIODICOFFSETS (
        RPBP_ESTIMATEMETAGENEBAYESFACTORS.out.bayes_factors
    )

    ch_extract_in = ch_bam
        .join(RPBP_SELECTPERIODICOFFSETS.out.periodic, by: 0)

    RPBP_EXTRACTORFPROFILES (
        ch_extract_in,
        ch_orfs_genomic_bed,
        ch_orfs_exons_bed
    )

    RPBP_ESTIMATEORFBAYESFACTORS (
        RPBP_EXTRACTORFPROFILES.out.profiles,
        ch_orfs_genomic_bed
    )

    RPBP_SELECTFINALPREDICTIONSET (
        RPBP_ESTIMATEORFBAYESFACTORS.out.bayes_factors,
        ch_genome_fasta
    )

    emit:
    metagene        = RPBP_EXTRACTMETAGENEPROFILES.out.metagene
    metagene_bf     = RPBP_ESTIMATEMETAGENEBAYESFACTORS.out.bayes_factors
    periodic        = RPBP_SELECTPERIODICOFFSETS.out.periodic
    orf_profiles    = RPBP_EXTRACTORFPROFILES.out.profiles
    lengths_offsets = RPBP_EXTRACTORFPROFILES.out.lengths_offsets
    orf_bayes       = RPBP_ESTIMATEORFBAYESFACTORS.out.bayes_factors
    predicted       = RPBP_SELECTFINALPREDICTIONSET.out.predicted
    dna_fasta       = RPBP_SELECTFINALPREDICTIONSET.out.dna_fasta
    protein_fasta   = RPBP_SELECTFINALPREDICTIONSET.out.protein_fasta
    versions        = ch_versions
}
