//
// Run riboWaltz P-site analysis and prepare MultiQC outputs
//

include { RIBOWALTZ     } from '../../../modules/nf-core/ribowaltz/main'
include { RIBOWALTZ_MQC } from '../../../modules/local/ribowaltz_mqc/main'

workflow RIBOWALTZ_QC {
    take:
    ch_transcriptome_bam    // channel: [ val(meta), path(bam) ]
    ch_gtf                  // channel: [ val(meta), path(gtf) ]
    ch_fasta                // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run riboWaltz for P-site offset calculation and QC
    //
    RIBOWALTZ(
        ch_transcriptome_bam,
        ch_gtf,
        ch_fasta
    )
    ch_versions = ch_versions.mix(RIBOWALTZ.out.versions)

    //
    // Prepare ribowaltz QC data for MultiQC custom content
    //
    RIBOWALTZ.out.ribowaltz_qc_data
        .map { meta, files -> files }
        .flatten()
        .branch { file ->
            psite_region: file.name.endsWith('.psite_region.tsv')
            frames: file.name.endsWith('.frames.tsv') && !file.name.contains('stratified')
            metaprofile: file.name.endsWith('.metaprofile_psite.tsv')
        }
        .set { ch_ribowaltz_qc }

    //
    // MODULE: Transform riboWaltz outputs to MultiQC custom content format
    //
    RIBOWALTZ_MQC(
        ch_ribowaltz_qc.psite_region.collect(),
        ch_ribowaltz_qc.frames.collect(),
        ch_ribowaltz_qc.metaprofile.collect()
    )

    ch_multiqc_files = ch_multiqc_files
        .mix(RIBOWALTZ_MQC.out.psite_regions_mqc.ifEmpty([]))
        .mix(RIBOWALTZ_MQC.out.frames_mqc.ifEmpty([]))
        .mix(RIBOWALTZ_MQC.out.metaprofile_start_mqc.ifEmpty([]))
        .mix(RIBOWALTZ_MQC.out.metaprofile_stop_mqc.ifEmpty([]))

    emit:
    offset          = RIBOWALTZ.out.offset              // channel: [ val(meta), path(offset) ]
    best_offset     = RIBOWALTZ.out.best_offset         // channel: [ val(meta), path(best_offset) ]
    psites          = RIBOWALTZ.out.psites              // channel: [ val(meta), path(psites) ]
    qc_plots        = RIBOWALTZ.out.ribowaltz_qc        // channel: [ val(meta), path(plots) ]
    qc_data         = RIBOWALTZ.out.ribowaltz_qc_data   // channel: [ val(meta), path(data) ]
    multiqc_files   = ch_multiqc_files                  // channel: [ path(mqc_files) ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}
