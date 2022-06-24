//
// Index the reference genome and align reads
//

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'

workflow BWAMEM2_ALIGNER {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    genome          // channel: file(ref_genome)

    main:

    ch_versions      = Channel.empty()
    ch_bwamem2_index = Channel.empty()

    //
    // index the genome with bwamem2 index
    //
    if (!params.bwamem2_index){
    BWAMEM2_INDEX ( genome )
    ch_bwamem2_index = BWAMEM2_INDEX.out.index
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else{
        ch_bwamem2_index = Channel.fromPath(params.bwamem2_index)
    }

    //
    // Map reads with bwamem2 mem
    //
    sort_bam = false
    BWAMEM2_MEM ( reads, ch_bwamem2_index, sort_bam )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    emit:
    bam            = BWAMEM2_MEM.out.bam  // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions          // channel: [ versions.yml ]
}
