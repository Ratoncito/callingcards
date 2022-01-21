//
// Check input samplesheet and get read channels
//

include { CREATE_PILEUP_DB } from '../../modules/local/create_pileup_db'
include { PROMOTER_ENRICHMENT } from '../../modules/local/promoter_enrichment'

workflow QUANTIFY_HOPS {
    take:
    pileup      // channel: [ val(meta), [ pileup ] ]
    barcode_length  // a value, length of the barcode, eg 13
    promoter_bed // path to the bedfile which describes the promoter regions
    background_data // path to the background data
    pileup_stranded // whether the enrichment should consider reads only on the same strand as the feature

    main:

    ch_versions = Channel.empty()

    //
    // create the pileup db
    //
    CREATE_PILEUP_DB ( pileup, barcode_length, background_data  )
    ch_versions = ch_versions.mix(CREATE_PILEUP_DB.out.versions)

    //
    // Calculate enrichment over promoter regions
    //
    PROMOTER_ENRICHMENT( CREATE_PILEUP_DB.out.pileup_db, promoter_bed, pileup_stranded  )
    ch_versions = ch_versions.mix(PROMOTER_ENRICHMENT.out.versions)

    emit:
    pileup_db = CREATE_PILEUP_DB.out.pileup_db                // channel: [ val(meta), path(pileup.sqlite) ]
    promoter_enrichment = PROMOTER_ENRICHMENT.out.enrichment  // channel: [ val(meta), path(promoter.csv) ]
    versions = ch_versions                                    // channel: [ versions.yml ]
}
