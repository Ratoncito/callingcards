
//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
// COPIED FROM nf-co/rnaseq
//

include { SAMTOOLS_SORT_INDEX_STATS } from '../nf-core/samtools_bam_sort_index_stats'
include { ADD_RG_AND_TAGS           } from '../../modules/local/add_read_group_and_tags'
include { COUNT_HOPS                } from '../../modules/local/count_hops'

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

    COUNT_HOPS(
        ADD_RG_AND_TAGS.out.bam,
        ADD_RG_AND_TAGS.out.bai,
        params.min_mapq
    )
    ch_version = ch_versions.mix(COUNT_HOPS.out.versions)


    emit:
    bam_index = ADD_RG_AND_TAGS.out.bam                // channel: [ val(meta), [ bam ], [ bai ] ]
    stats     = SAMTOOLS_SORT_INDEX_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat  = SAMTOOLS_SORT_INDEX_STATS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats  = SAMTOOLS_SORT_INDEX_STATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    bed       = COUNT_HOPS.out.bed                     // channel: [ val(meta), [ bed ] ]
    versions  = ch_versions                            // channel: [ versions.yml ]
}
