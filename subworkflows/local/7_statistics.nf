//
// Check input samplesheet and get read channels
//

include { YEAST_FIND_SIG_PROMOTERS } from "${projectDir}/modules/local/yeast/find_sig_promoters"

workflow STATISTICS {
    take:
    qbed      // channel: [ val(meta), file(bed) ]
    promoter_bed
    background_qbed
    chr_map
    standard_chr_format
    sqlite_db_out
    poisson_pseudocount

    main:

    ch_versions             = Channel.empty()
    ch_sig_promoters        = Channel.empty()
    ch_sig_promoters_sqlite = Channel.empty()

    if(params.organism == 'yeast'){
        YEAST_FIND_SIG_PROMOTERS (
            qbed,
            promoter_bed,
            background_qbed,
            chr_map,
            standard_chr_format,
            sqlite_db_out,
            poisson_pseudocount
        )
        ch_sig_promoters.mix(YEAST_FIND_SIG_PROMOTERS.out.sig_promoters)
        ch_sig_promoters_sqlite.mix(YEAST_FIND_SIG_PROMOTERS.out.sig_promoters_sqlite)
        ch_versions = ch_versions.mix(YEAST_FIND_SIG_PROMOTERS.out.versions)
    } else {
        println 'Statistics module not implemented for your organism'
    }

    emit:
    sig_promoters        = ch_sig_promoters                 // channel: [ val(meta), csv ]
    sig_promoters_sqlite = ch_sig_promoters_sqlite          // channel: [ val(meta), sqlite_db ]
    versions             = ch_versions                      // channel: [ versions.yml ]
}
