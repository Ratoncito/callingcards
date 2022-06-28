//
// Check input samplesheet and get read channels
//

include { BARCODE_QC } from '../../modules/local/barcode_qc'

workflow PROCESS_HOPS {
    take:
    bed      // channel: [ val(meta), file(bed) ]
    barcode_details // channel: [val(meta), file(barcode_details)]

    main:

    ch_versions = Channel.empty()

    //
    // parse the 'name' column of the bed 6 x ? calling cards bed, calculate
    // position probability matricies and tallys of the varieties of seqeunces
    // filter rows of the bed file down to those rows which meet expectations on
    // barcode (and insert sequence). If a tf_map is provided, split the bed
    // file into individual TFs
    //
    BARCODE_QC ( bed, barcode_details )
    ch_versions = ch_versions.mix(BARCODE_QC.out.versions)

    emit:
    bed   = BARCODE_QC.out.bed    // channel: [ val(meta), [bed file(s)] ]
    ppm   = BARCODE_QC.out.ppm    // channel: [ file(position prob matrix.tsv) ]
    tally = BARCODE_QC.out.tally  // channel: [ file(tally.tsv) ]
    versions = ch_versions        // channel: [ versions.yml ]
}
