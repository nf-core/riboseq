//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_CANONICAL_GTF    } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { GFFREAD as GFFREAD_CANONICAL      } from '../../../modules/nf-core/gffread'
include { AGAT_SPKEEPLONGESTISOFORM         } from '../../../modules/nf-core/agat/spkeeplongestisoform'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { SALMON_INDEX as SALMON_INDEX_TE   } from '../../../modules/nf-core/salmon/index'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    salmon_index             // directory: /path/to/salmon/index/
    bbsplit_index            // directory: /path/to/rsem/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star'
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             //   boolean: Skip BBSplit for removal of non-reference genome reads
    skip_sortmerna           //   boolean: Skip sortmerna for removal of non-reference genome reads
    ribo_removal_tool        //    string: Tool for rRNA removal ('sortmerna', 'bowtie2', or 'ribodetector')
    skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline
    build_te_pseudo_index    //   boolean: Build Salmon index for TE pseudo-alignment
    canonical_gtf            //      file: /path/to/canonical.gtf (one-transcript-per-gene backbone; null to derive)

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
            } else {
                ch_gtf = Channel.value(file(gtf))
            }
        } else if (gff) {
            def ch_gff
            if (gff.endsWith('.gz')) {
                ch_gff = GUNZIP_GFF ( [ [:], gff ] ).gunzip
            } else {
                ch_gff = Channel.value(file(gff)).map { item -> [ [:], item ] }
            }
            ch_gtf = GFFREAD ( ch_gff, [] ).gtf.map { it[1] }
        }

        // Determine whether to filter the GTF or not
        def filter_gtf =
            ((
                // Condition 1: Alignment is required and aligner is set
                !skip_alignment && aligner
            ) ||
            (
                // Condition 2: Transcript FASTA file is not provided
                !transcript_fasta
            )) &&
            (
                // Condition 4: --skip_gtf_filter is not provided
                !skip_gtf_filter
            )
        if (filter_gtf) {
            GTF_FILTER ( ch_fasta, ch_gtf )
            ch_gtf = GTF_FILTER.out.genome_gtf
            ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        }
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : "gene_biotype"
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], additional_fasta ] ).gunzip.map { it[1] }
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta))
        }

        CUSTOM_CATADDITIONALFASTA(
            ch_fasta.combine(ch_gtf).map{ fa, gt -> [[:], fa, gt] },
            ch_add_fasta.map{[[:], it]},
            biotype
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map{it[1]}.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map{it[1]}.first()
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    //
    // Build the canonical (one-transcript-per-gene) annotation backbone used by
    // ORF callers, riboWaltz, plastid P-site quantification and DTE. The full
    // multi-isoform `ch_gtf` is still used for genome-guided alignment.
    //
    ch_canonical_gtf = ch_gtf
    if (canonical_gtf) {
        if (canonical_gtf.endsWith('.gz')) {
            ch_canonical_gtf = GUNZIP_CANONICAL_GTF ( [ [:], canonical_gtf ] ).gunzip.map { it[1] }
        } else {
            ch_canonical_gtf = Channel.value(file(canonical_gtf))
        }
    } else if (gtf || gff) {
        // No explicit canonical GTF: derive longest-isoform per gene with AGAT.
        // Versions are emitted via the `versions` topic channel and collected
        // workflow-wide; no explicit mix into ch_versions required.
        AGAT_SPKEEPLONGESTISOFORM (
            ch_gtf.map { [ [id: it.baseName], it ] },
            []
        )
        GFFREAD_CANONICAL (
            AGAT_SPKEEPLONGESTISOFORM.out.gff,
            []
        )
        ch_canonical_gtf = GFFREAD_CANONICAL.out.gtf.map { it[1] }
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
    }

    //
    // Create chromosome sizes file
    //
    SAMTOOLS_FAIDX ( ch_fasta.map { [ [:], it, [] ] }, true )
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { it[1] }
    //
    // Get list of indices that need to be created
    //
    def prepare_tool_indices = []
    if (!skip_bbsplit) { prepare_tool_indices << 'bbsplit' }
    if (!skip_sortmerna) { prepare_tool_indices << 'sortmerna' }
    if (!skip_alignment) { prepare_tool_indices << aligner }

    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ [:], bbsplit_index ] ).untar.map { it[1] }
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index))
            }
        } else {
            Channel
                .from(file(bbsplit_fasta_list))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fa -> [ [ 'id', id ], [ 'fasta', file(fa, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { it -> it[1] } // Get rid of keys and keep grouped values
                .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT ( [ [:], [] ], [], ch_fasta, ch_bbsplit_fasta_list, true ).index
        }
    }

    //
    // Prepare rRNA fastas for rRNA removal (sortmerna, bowtie2, or ribodetector)
    //
    ch_sortmerna_index = Channel.empty()
    ch_rrna_fastas = Channel.empty()

    // Populate ch_rrna_fastas when sortmerna or bowtie2 is selected (ribodetector uses its own model)
    if (ribo_removal_tool in ['sortmerna', 'bowtie2']) {
        ribo_db = file(sortmerna_fasta_list)
        ch_rrna_fastas = Channel.from(ribo_db.readLines())
            .map { row -> file(row, checkIfExists: true) }
    }

    // Only build sortmerna index when sortmerna is selected as the rRNA removal tool
    if ('sortmerna' in prepare_tool_indices) {
        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ( [ [:], sortmerna_index ] ).untar.map { it[1] }
            } else {
                ch_sortmerna_index = Channel.value(file(sortmerna_index))
            }
        } else {
            SORTMERNA_INDEX (
                Channel.of([ [],[] ]),
                ch_rrna_fastas.collect().map { [ 'rrna_refs', it ] },
                Channel.of([ [],[] ])
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
        } else {
            // Check if an AWS iGenome has been provided to use the appropriate version of STAR
            def is_aws_igenome = false
            if (fasta && gtf) {
                if ((file(fasta).getName() - '.gz' == 'genome.fa') && (file(gtf).getName() - '.gz' == 'genes.gtf')) {
                    is_aws_igenome = true
                }
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
            }
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
        }
    }

    //
    // Build Salmon index for TE pseudo-alignment (with lower k-mer size for short reads)
    //
    ch_salmon_index_te = Channel.empty()
    if (build_te_pseudo_index) {
        ch_salmon_index_te = SALMON_INDEX_TE ( ch_fasta, ch_transcript_fasta ).index
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    canonical_gtf    = ch_canonical_gtf          // channel: path(canonical.gtf) - one-transcript-per-gene backbone
    fai              = ch_fai                    // channel: path(genome.fai)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    rrna_fastas      = ch_rrna_fastas            // channel: path(sortmerna_fasta_list)
    sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    salmon_index_te  = ch_salmon_index_te        // channel: path(salmon_te/index/) - for TE pseudo-alignment
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
