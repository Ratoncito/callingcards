//
// Check input samplesheet and get read channels
//

include { BARCODE_QC_DEMULTIPLEX } from '../../modules/local/barcodeQC_demultiplex'
include { YEAST_FIND_SIG_PROMOTERS } from '../../modules/local/yeast_find_sig_promoters'

workflow PROCESS_HOPS {
    take:
    bed      // channel: [ val(meta), file(bed) ]
    barcode_details // channel: [val(meta), file(barcode_details)]

    main:

    ch_versions = Channel.empty()
    ch_sig_promoters = Channel.empty()

    //
    // parse the 'name' column of the bed 6 x ? calling cards bed, calculate
    // position probability matricies and tallys of the varieties of seqeunces
    // filter rows of the bed file down to those rows which meet expectations on
    // barcode (and insert sequence). If a tf_map is provided, split the bed
    // file into individual TFs
    //
    BARCODE_QC_DEMULTIPLEX ( bed, barcode_details )
    ch_versions = ch_versions.mix(BARCODE_QC_DEMULTIPLEX.out.versions)

    if(params.organism == 'yeast'){
        YEAST_FIND_SIG_PROMOTERS (
            BARCODE_QC_DEMULTIPLEX.out.whatever
        )
        ch_sig_promoters.mix(YEAST_FIND_SIG_PROMOTERS.out.sig_promoters)
    } else if(params.organism == 'mammal'){

    } else{
        exit 1, 'Find Significant Promoters: Organism not recognized'
    }

    emit:
    //sig_promoters   = ch_sig_promoters    // channel: [ val(meta), [bed file(s)] ]
   // ppm   = BARCODE_QC_DEMULTIPLEX.out.ppm    // channel: [ file(position prob matrix.tsv) ]
    //tally = BARCODE_QC_DEMULTIPLEX.out.tally  // channel: [ file(tally.tsv) ]
    versions = ch_versions        // channel: [ versions.yml ]
}
