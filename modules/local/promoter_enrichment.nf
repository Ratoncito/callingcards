process PROMOTER_ENRICHMENT {

    tag "$meta.id"
    label 'process_medium'

    container "library://cmatkhan/default/calling_cards_tools:sha256.6987abea67e7fcef44443fff964de003be559f56c839ebf43c2d4160a92b675a"

    input:
    tuple val(meta), path(pileup_db)
    path  promoter_bed
    val   stranded

    output:
    tuple val(meta), path("*.csv"), emit: enrichment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/callingcards/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    promoter_enrichment.R \\
        -p ${pileup_db} \\
        -b ${promoter_bed} \\
        -m ${meta.barcodes} \\
        -s ${stranded} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
