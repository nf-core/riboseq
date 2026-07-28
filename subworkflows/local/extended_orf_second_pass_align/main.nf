//
// Extended-ORF second STAR pass against the hybrid transcriptome.
//
// RiboCode requires a transcriptome-coordinate BAM keyed to whichever
// transcriptome FASTA was used at alignment time. To bring novel intergenic
// transcripts into RiboCode, this subworkflow:
//   1. Rebuilds the transcriptome FASTA from the hybrid GTF (gffread -w).
//   2. Rebuilds a STAR index against the genome FASTA + hybrid GTF.
//   3. Aligns Ribo-seq reads only against the hybrid STAR index, emitting a
//      hybrid transcriptome BAM that RiboCode can consume.
//
// Compute cost: roughly doubles STAR alignment work for Ribo-seq samples.
// Gating (params.extended_orf_analysis + a configured novel source) lives at
// the call site.
//

include { GFFREAD as GFFREAD_HYBRID_TRANSCRIPTOME           } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HYBRID } from '../../../modules/nf-core/star/genomegenerate'
include { FASTQ_ALIGN_STAR as FASTQ_ALIGN_STAR_HYBRID       } from '../../nf-core/fastq_align_star'

workflow EXTENDED_ORF_SECOND_PASS_ALIGN {

    take:
    // The two reference channels below must be value channels: they are paired with the per-sample
    // read channel in the STAR fan-out, and a single-item queue would serve only one sample.
    ch_hybrid_gtf            // value channel: path(hybrid.gtf)
    ch_fasta                 // value channel: path(genome.fasta)
    ch_reads_for_alignment   // channel: [ val(meta), [ path(fastq), ... ] ]
    star_ignore_sjdbgtf      // val: params.star_ignore_sjdbgtf

    main:

    // Extract spliced transcript sequences (canonical + novel) from the hybrid
    // GTF, then rebuild a STAR index against the same genome FASTA + hybrid
    // GTF so a second STAR pass can emit a hybrid transcriptome BAM.
    GFFREAD_HYBRID_TRANSCRIPTOME(
        ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] },
        ch_fasta
    )

    STAR_GENOMEGENERATE_HYBRID(
        ch_fasta.map { fasta -> [ [id: 'hybrid_reference'], fasta ] },
        ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] }
    )

    ch_hybrid_transcriptome_fasta = GFFREAD_HYBRID_TRANSCRIPTOME.out.gffread_fasta.map { _meta, fasta -> fasta }
    ch_hybrid_star_index          = STAR_GENOMEGENERATE_HYBRID.out.index.map { _meta, index -> index }

    // Only Ribo-seq samples feed RiboCode; skip the costly re-alignment for
    // RNA-seq and TI-seq.
    ch_reads_for_hybrid_alignment = ch_reads_for_alignment
        .filter { meta, _reads -> meta.sample_type == 'riboseq' }

    FASTQ_ALIGN_STAR_HYBRID(
        ch_reads_for_hybrid_alignment,
        ch_hybrid_star_index.map { [ [:], it ] },
        ch_hybrid_gtf.map { [ [:], it ] },
        star_ignore_sjdbgtf,
        ch_fasta.map { [ [:], it, [] ] },
        ch_hybrid_transcriptome_fasta.map { [ [:], it, [] ] }
    )

    ch_multiqc_files = FASTQ_ALIGN_STAR_HYBRID.out.stats.collect{it[1]}
        .mix(FASTQ_ALIGN_STAR_HYBRID.out.flagstat.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR_HYBRID.out.idxstats.collect{it[1]})
        .mix(FASTQ_ALIGN_STAR_HYBRID.out.log_final.collect{it[1]})

    emit:
    transcriptome_bam = FASTQ_ALIGN_STAR_HYBRID.out.orig_bam_transcript // [ meta, bam ]
    genome_bam        = FASTQ_ALIGN_STAR_HYBRID.out.bam                  // [ meta, bam ] - coordinate-sorted
    genome_bai        = FASTQ_ALIGN_STAR_HYBRID.out.index                // [ meta, bai ]
    transcript_fasta  = ch_hybrid_transcriptome_fasta                   // channel: path(fasta)
    multiqc_files     = ch_multiqc_files                                // channel: path(file)
}
