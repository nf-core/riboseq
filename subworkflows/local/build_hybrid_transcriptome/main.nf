//
// Build a hybrid transcriptome FASTA + STAR index for the extended-ORF
// second-pass alignment (issue #171).
//
// STAR `--quantMode TranscriptomeSAM` produces a transcriptome-coordinate
// BAM keyed to the transcriptome FASTA supplied at alignment time. Novel
// StringTie transcripts do not exist in the reference transcriptome FASTA
// and cannot be counted from the primary transcriptome BAM. This subworkflow
// extracts spliced transcript sequences from the hybrid GTF (canonical +
// filtered novel intergenic) and rebuilds a STAR index against the same
// genome FASTA + hybrid GTF so a second STAR pass can emit a hybrid
// transcriptome BAM that RiboCode can consume.
//
// riboWaltz is intentionally NOT switched onto the hybrid transcriptome.
// riboWaltz is a QC/calibration tool; its frame plots and metaheatmaps are
// driven by CDS-bearing canonical transcripts. Feeding it CDS-absent novel
// transcripts would degrade diagnostic plots without any ORF-discovery gain.
// See docs/usage.md (Extended ORF discovery) for the full rationale.
//

include { GFFREAD as GFFREAD_HYBRID_TRANSCRIPTOME } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HYBRID } from '../../../modules/nf-core/star/genomegenerate'

workflow BUILD_HYBRID_TRANSCRIPTOME {

    take:
    ch_fasta        // channel: path(genome.fasta)
    ch_hybrid_gtf   // channel: path(hybrid_reference.gtf)

    main:

    //
    // Extract spliced transcript sequences (canonical + novel) from the
    // hybrid GTF into a single transcriptome FASTA. `-w` writes the spliced
    // exon FASTA; gffread requires the genome FASTA via `-g`.
    //
    GFFREAD_HYBRID_TRANSCRIPTOME(
        ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] },
        ch_fasta
    )

    //
    // Rebuild the STAR index against the same genome FASTA but with the
    // hybrid GTF as `--sjdbGTFfile`. STAR uses the GTF only to seed splice
    // junctions and the transcriptome assignment table; the genome SA index
    // itself is identical to the primary index, so re-aligning against this
    // index is what gets us hybrid transcriptome coordinates in the resulting
    // `Aligned.toTranscriptome.out.bam`.
    //
    STAR_GENOMEGENERATE_HYBRID(
        ch_fasta.map { fasta -> [ [id: 'hybrid_reference'], fasta ] },
        ch_hybrid_gtf.map { gtf -> [ [id: 'hybrid_reference'], gtf ] }
    )

    emit:
    hybrid_transcriptome_fasta = GFFREAD_HYBRID_TRANSCRIPTOME.out.gffread_fasta.map { _meta, fasta -> fasta }
    hybrid_star_index          = STAR_GENOMEGENERATE_HYBRID.out.index.map { _meta, index -> index }
}
