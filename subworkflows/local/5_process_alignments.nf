
//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
// COPIED FROM nf-co/rnaseq
//

include { SAMTOOLS_SORT_INDEX_STATS } from '../nf-core/samtools_bam_sort_index_stats'
include { ADD_RG_AND_TAGS           } from '../../modules/local/add_read_group_and_tags'
include { SAMTOOLS_MPILEUP          } from '../../modules/nf-core/modules/samtools/mpileup/main'

workflow PROCESS_ALIGNMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    fasta // path(genome.fasta) path to the fasta file
    fai // channel: [val(meta), path(fasta index)] note that the meta is empty

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_SORT_INDEX_STATS(
        bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_STATS.out.versions)

    ADD_RG_AND_TAGS (
        SAMTOOLS_SORT_INDEX_STATS.out.bam_index,
        fasta,
        fai,
        params.barcode_length,
        params.insertion_length
    )
    ch_version = ch_versions.mix(ADD_RG_AND_TAGS.out.versions)

    SAMTOOLS_MPILEUP(
        ADD_RG_AND_TAGS.out.bam_index,
        fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)


    emit:
    bam_index = ADD_RG_AND_TAGS.out.bam_index // channel: [ val(meta), [ bam ], [ bai ] ]
    stats     = SAMTOOLS_SORT_INDEX_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat  = SAMTOOLS_SORT_INDEX_STATS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats  = SAMTOOLS_SORT_INDEX_STATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    mpileup   = SAMTOOLS_MPILEUP.out.mpileup           // channel: [ val(meta), [ mpileup ] ]
    versions  = ch_versions                   // channel: [ versions.yml ]

}
