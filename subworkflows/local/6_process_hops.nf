//
// Check input samplesheet and get read channels
//

include { BARCODE_QC_DEMULTIPLEX } from "${projectDir}/modules/local/barcodeQC_demultiplex"
include { YEAST_FIND_SIG_PROMOTERS } from "${projectDir}/modules/local/yeast/find_sig_promoters"

workflow PROCESS_HOPS {
    take:
    bed      // channel: [ val(meta), file(bed) ]
    barcode_details // channel: [val(meta), file(barcode_details)]
    promoter_bed
    background_ccf
    chr_map
    standard_chr_format
    sqlite_db_out
    poisson_pseudocount

    main:

    ch_versions             = Channel.empty()
    ch_stats                = Channel.empty()
    ch_sig_promoters        = Channel.empty()
    ch_sig_promoters_sqlite = Channel.empty()

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
            BARCODE_QC_DEMULTIPLEX.out.ccf,
            promoter_bed,
            background_ccf,
            chr_map,
            standard_chr_format,
            sqlite_db_out,
            poisson_pseudocount
        )
        ch_sig_promoters.mix(YEAST_FIND_SIG_PROMOTERS.out.sig_promoters)
        ch_sig_promoters_sqlite.mix(YEAST_FIND_SIG_PROMOTERS.out.sig_promoters_sqlite)
    } else if(params.organism == 'mammal'){
        exit 1, 'Process Hops Error: NotImplementedError: mammals'
    } else{
        exit 1, 'Process Hops Error:  Organism not recognized'
    }

    emit:
    sig_promoters = ch_sig_promoters                 // channel: [ val(meta), csv ]
    sig_promoters_sqlite = ch_sig_promoters_sqlite   // channel: [ val(meta), sqlite_db ]
    ppm           = BARCODE_QC_DEMULTIPLEX.out.ppm   // channel: [ file(position prob matrix.tsv) ]
    tally         = BARCODE_QC_DEMULTIPLEX.out.tally // channel: [ file(tally.tsv) ]
    versions      = ch_versions                      // channel: [ versions.yml ]
}
