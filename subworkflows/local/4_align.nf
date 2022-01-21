//
// Align reads to a reference genome
// note that this can be parameterized -- could put $param.aligner
// in the include ... from ... path below
//

include { BWAMEM2_ALIGNER } from '../nf-core/bwamem2'

workflow ALIGN {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    genome           // channel: file(ref_genome)

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    if(params.aligner == 'bwamem2') {
        BWAMEM2_ALIGNER (
            reads,
            genome
        )
        ch_bam      = ch_bam.mix(BWAMEM2_ALIGNER.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_ALIGNER.out.versions)
    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }


    emit:
    bam       = ch_bam           // channel: [ val(meta), [ bam ] ]
    versions  = ch_versions      // channel: [ versions.yml ]
}
