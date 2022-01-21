process CREATE_PILEUP_DB {

    tag "$meta.id"
    label 'process_medium'

    container "library://cmatkhan/default/calling_cards_tools:sha256.6987abea67e7fcef44443fff964de003be559f56c839ebf43c2d4160a92b675a"

    input:
    tuple val(meta), path(pileup)
    val   barcode_length
    path  background_data

    output:
    tuple val(meta), path("*.sqlite"), emit: pileup_db
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_pileup_db.R -p ${pileup} \\
                       -b ${barcode_length} \\
                       -o ${prefix} \\
                       -d ${background_data}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
