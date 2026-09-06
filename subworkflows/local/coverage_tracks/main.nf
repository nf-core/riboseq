//
// Strand-aware genome-coverage tracks for IGV/JBrowse.
//
// For stranded libraries the genome BAMs are first split into per-strand BAMs
// via a samtools view flag filter (forward / reverse). Unstranded libraries
// are passed straight through. Each resulting BAM becomes a bedgraph
// (bedtools genomecov), then a bigWig (UCSC bedGraphToBigWig).
//

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SPLIT_BY_STRAND } from '../../../modules/nf-core/samtools/view'
include { BEDTOOLS_GENOMECOV                             } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG                          } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

// A filter for samtools view which splits alignments by first-of-pair strand,
// taking into consideration the strandedness of the library. Used by
// SAMTOOLS_VIEW_SPLIT_BY_STRAND.
def getStrandFilter(strandedness, strand) {
    def sameOrientation = (strand == 'forward') == (strandedness == 'forward')
    sameOrientation
        ? "-e '((flag.read1 || !flag.paired) && !flag.reverse) || (flag.read2 &&  flag.reverse)'"
        : "-e '((flag.read1 || !flag.paired) &&  flag.reverse) || (flag.read2 && !flag.reverse)'"
}

workflow COVERAGE_TRACKS {

    take:
    ch_genome_bam       // channel: [ val(meta), path(bam) ]
    ch_genome_bam_index // channel: [ val(meta), path(bai) ]
    ch_fai              // channel: path(genome.fai)

    main:

    // Stranded libraries: split BAMs by mate-1 strand before coverage.
    ch_split_by_strand = ch_genome_bam
        .join(ch_genome_bam_index, by: [0])
        .filter { meta, _bam, _bai -> meta.strandedness in ['forward', 'reverse'] }

    SAMTOOLS_VIEW_SPLIT_BY_STRAND(
        ch_split_by_strand
            .flatMap { meta, bam, bai ->
                ['forward', 'reverse'].collect { strand ->
                    def strand_filter = getStrandFilter(meta.strandedness, strand)
                    [meta + [strand: strand, strand_filter: strand_filter], bam, bai]
                }
            },
        [[], [], []],  // No reference fasta/fai
        [[], []],      // No qname file
        [[], []],      // No region BED
        []             // No index format
    )

    // bedgraphs from the per-strand BAMs and from any unstranded BAMs.
    BEDTOOLS_GENOMECOV(
        SAMTOOLS_VIEW_SPLIT_BY_STRAND.out.bam
            .map { meta, bam -> [meta, bam, 1] }
            .mix(ch_genome_bam
                .filter { meta, _bam -> meta.strandedness == 'unstranded' }
                .map { meta, bam -> [meta + [strand: 'unstranded'], bam, 1] }
            ),
        ch_fai,
        'bedgraph',
        false
    )

    UCSC_BEDGRAPHTOBIGWIG(
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_fai
    )

    emit:
    bigwig = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), path(bigwig) ]
}
