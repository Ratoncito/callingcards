process YEAST_FIND_SIG_PROMOTERS {

    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.9.13 conda-forge::pandas=1.4.3 conda-forge::scipy=1.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a1289c2d7470e63e3c3a9f6131984bbf7c28ad45:92387cc4b89149a92d9f98955d1624b48f915ddb-0' :
        'quay.io/biocontainers/mulled-v2-a1289c2d7470e63e3c3a9f6131984bbf7c28ad45:92387cc4b89149a92d9f98955d1624b48f915ddb-0' }"

    input:
    tuple val(meta), path(expr_ccf)
    path(promoter_bed)
    path(background_ccf)
    path(chr_map)
    val standard_chr_format
    val sqlite_db_out
    val poisson_pseudocount


    output:
    tuple val(meta), path("*_stats.csv"), emit: sig_promoters // [val(meta), path(csv) ]
    tuple val(meta), path("*.sqlite"), optional: true,    emit: sig_promoters_sqlite // [ val(meta), path(sqlite_db)]
    path  "versions.yml",                 emit: versions

    script: // see nf-core-callingcards/bin/mammals_barcode_qc.py
    """
        yeast_find_sig_promoters.py \
            -p $promoter_bed \
            -b $background_ccf \
            -e $expr_ccf \
            -c $chr_map \
            -s $standard_chr_format \
            -d $sqlite_db_out \
            -x $poisson_pseudocount

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip freeze | grep pandas | sed 's/pandas==//g')
        scipy: \$(pip freeze | grep scipy | sed 's/scipy==//g')
    END_VERSIONS

    """
}
