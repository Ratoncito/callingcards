//
// Check input samplesheet and get read channels
//

// TODO handle barcode details differently?
// TODO check that sample ids are unique

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .multiMap { it ->
            reads: create_fastq_channel(it)
            barcode_details: create_barcode_details_channel(it)}
        .set { out }

    emit:
    reads = out.reads                         // channel: [ val(meta), [ file(r1 fastq), file(f2_fastq)* depending on input ]]
    barcode_details = out.barcode_details     // channel: [ val(meta), file(barcode_details.json) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
// where meta contains the sample id and path to the barcodes file
// dummy_file is a string to a file which exists, but is empty
// cite: nf-core/rnaseq pipeline
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    def fastq_1 = file(row.fastq_1).exists() ?
                    row.fastq_1 :
                    "${projectDir}/assets/test_data/${row.fastq_1}"

    if (!file(fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(fastq_1) ]]
    } else {
        def fastq_2 = file(row.fastq_2).exists() ?
                        row.fastq_2 :
                        "${projectDir}/assets/test_data/${row.fastq_2}"

        if (!file(fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(fastq_1), file(fastq_2) ] ]
    }
    return fastq_meta
}

// Like create_fastq_channel, this parses a row of the input sample sheet
// and returns an array where the first item is a map of metadata,
// and the second item is the barcode_details.json as a file object
def create_barcode_details_channel(LinkedHashMap row){
    def meta = [:]
    def barcode_details_arr = []
    meta.id = row.sample

    def barcode_details = file(row.barcode_details).exists() ?
                    row.barcode_details :
                    "${projectDir}/assets/${row.barcode_details}"

    if (!file(barcode_details).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Barcode detail does not exist!\n${row.barcode_details}"
    }

    barcode_details_arr = [meta, [file(barcode_details)]]

    return barcode_details_arr
}
