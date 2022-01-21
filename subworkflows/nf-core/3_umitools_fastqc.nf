//
// Read QC, UMI extraction and trimming
// NOTE: copied from nf-co/rnaseq/subworkflows/nf-core/fastqc_umitools_trimgalore.nf
//

// TODO this is the longest process -- split into n groups and do in parallel
// use this as template https://www.nextflow.io/example3.html

include { FASTQC           } from '../../modules/nf-core/modules/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/modules/umitools/extract/main'

workflow UMITOOLS_FASTQC {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    // todo error checking on this -- should not be empty?
    umi_log   = Channel.empty()
    umi_reads = Channel.empty()
    // run umi extract to add barcodes to fastq id lines
    UMITOOLS_EXTRACT ( reads ).reads.set { umi_reads }

    umi_log     = UMITOOLS_EXTRACT.out.log
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

    // run fastqc after trimming off the barcodes, etc
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()

    FASTQC ( reads ).html.set { fastqc_html }
    fastqc_zip  = FASTQC.out.zip
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:
    reads = umi_reads  // TODO channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
