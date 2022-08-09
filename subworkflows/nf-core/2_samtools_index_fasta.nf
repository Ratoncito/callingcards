//
// Index the reference genome and align reads
//

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'

workflow SAMTOOLS_INDEX_FASTA {
    take:
    fasta   // path(genome)

    main:

    ch_versions = Channel.empty()

    //
    // index the genome with bwamem2 index
    //
    ch_genome = Channel.of( [''] ).mix(fasta)
    SAMTOOLS_FAIDX ( ch_genome )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    fai      = SAMTOOLS_FAIDX.out.fai  // channel: path(*.fai)
    versions = ch_versions             // channel: [ versions.yml ]
}
