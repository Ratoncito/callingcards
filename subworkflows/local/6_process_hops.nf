//
// Process hop bed file to qbed by TF.
// also output a number of QC metrics
//

include { BARCODE_QC_DEMULTIPLEX } from "${projectDir}/modules/local/barcodeQC_demultiplex"

workflow PROCESS_HOPS {
    take:
    process_hops_input // [meta, file(barcode_details), [file(passing.bed), file(failing.bed)]]

    main:

    ch_versions             = Channel.empty()

    //
    // parse the 'name' column of the bed 6 x ? calling cards bed, calculate
    // position probability matricies and tallys of the varieties of seqeunces
    // filter rows of the bed file down to those rows which meet expectations on
    // barcode (and insert sequence). If a tf_map is provided, split the bed
    // file into individual TFs
    BARCODE_QC_DEMULTIPLEX ( process_hops_input )
    ch_versions = ch_versions.mix(BARCODE_QC_DEMULTIPLEX.out.versions)

    BARCODE_QC_DEMULTIPLEX.out.qbed
        .transpose()
        .map{ meta, qbed -> augment_meta(meta,qbed)}
        .set{ ch_demult_qbed }

    emit:

    qbed   = ch_demult_qbed                   // channel: [ val(meta), csv ]
    ppm    = BARCODE_QC_DEMULTIPLEX.out.ppm   // channel: [ file(position prob matrix.tsv) ]
    tally  = BARCODE_QC_DEMULTIPLEX.out.tally // channel: [ file(tally.tsv) ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

def augment_meta(Map meta, qbed) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta["sample_id"] = meta.id

    new_meta.id = qbed.baseName-"_bc_fltr"

    return [new_meta, file(qbed)]
}

