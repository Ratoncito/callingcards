//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
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



    meta.barcodes     = (row.barcodes == "TEST_INPUT.txt") ?
                           "${projectDir}/assets/mammals_barcode_safelist.txt" :
                           file(row.barcodes)

    if(!file(meta.barcodes).exists()){
        exit 1, "ERROR: Please check input samplesheet -> the barcodes file does not exist!\n${meta.barcodes}"

    }
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (row.fastq_1 != "TEST_INPUT.fq.gz" & !file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = (row.fastq_1 == "TEST_INPUT.fq.gz") ?
                        [ meta, file("${projectDir}/assets/test_data/AY53-1_50k_downsampled_human.fastq.gz") ] :
                         [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
