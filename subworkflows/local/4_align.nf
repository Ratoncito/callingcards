//
// Align reads to a reference genome
// note that this can be parameterized -- could put $param.aligner
// in the include ... from ... path below
//

include { BWAMEM2_MEM } from "${projectDir}/modules/nf-core/modules/bwamem2/mem/main.nf"

workflow ALIGN {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    bwamem2_index // channel: file(fasta)

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    if(params.aligner == 'bwamem2') {
        sort_bam = false
        BWAMEM2_MEM (
            reads,
            bwamem2_index,
            sort_bam
        )
        ch_bam      = ch_bam.mix(BWAMEM2_MEM.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }


    emit:
    bam       = ch_bam           // channel: [ val(meta), [ bam ] ]
    versions  = ch_versions      // channel: [ versions.yml ]
}
