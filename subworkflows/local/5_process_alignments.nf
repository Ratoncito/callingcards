
//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
// COPIED FROM nf-co/rnaseq
//

include { SAMTOOLS_SORT_INDEX_STATS     } from "${projectDir}/subworkflows/nf-core/samtools_bam_sort_index_stats"
include { ADD_RG_AND_TAGS               } from "${projectDir}/modules/local/add_read_group_and_tags"
include { PICARD_COLLECTMULTIPLEMETRICS } from "${projectDir}/modules/nf-core/modules/picard/collectmultiplemetrics/main"
include { PRESEQ_LCEXTRAP               } from "${projectDir}/modules/nf-core/modules/preseq/lcextrap/main"
include { PRESEQ_CCURVE                 } from "${projectDir}/modules/nf-core/modules/preseq/ccurve/main"
include { COUNT_HOPS                    } from "${projectDir}/modules/local/count_hops"

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

    PICARD_COLLECTMULTIPLEMETRICS(
        bam,
        fasta,
        fai.map{meta,fai -> fai}
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    PRESEQ_LCEXTRAP(
        bam
    )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)

    PRESEQ_CCURVE(
        bam
    )
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions)

    COUNT_HOPS(
        ADD_RG_AND_TAGS.out.bam,
        ADD_RG_AND_TAGS.out.bai,
        params.min_mapq
    )
    ch_version = ch_versions.mix(COUNT_HOPS.out.versions)


    emit:
    bam_index       = ADD_RG_AND_TAGS.out.bam                   // channel: [ val(meta), [ bam ], [ bai ] ]
    stats           = SAMTOOLS_SORT_INDEX_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat        = SAMTOOLS_SORT_INDEX_STATS.out.flagstat    // channel: [ val(meta), [ flagstat ] ]
    idxstats        = SAMTOOLS_SORT_INDEX_STATS.out.idxstats    // channel: [ val(meta), [ idxstats ] ]
    picard_metrics  = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ metrics ] ]
    preseq_lcextrap = PRESEQ_LCEXTRAP.out.lc_extrap             // channel: [ val(meta), [ metrics ] ]
    preseq_ccurve   = PRESEQ_CCURVE.out.c_curve                 // channel: [ val(meta), [ metrics ] ]
    bed             = COUNT_HOPS.out.bed                         // channel: [ val(meta), [ bed ] ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}
