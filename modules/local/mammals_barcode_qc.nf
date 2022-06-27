process MAMMAL_BARCODE_QC {

    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::pysam=0.17.0 pandas=1.4.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'quay.io/biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(bed)
    path(barcode_components)
    path(barcode_component_indicies)
    val insert_seq


    output:
    tuple val(meta), path("*_bc_fltr.bed"), emit: // [val(meta), path(bed) ]
    path("*_tally.tsv"),                  , emit: // [ val(mega), path(tally)]
    path("*_ppm.tsv"),                    , emit: // [ val(mega), path(ppm)]
    path  "versions.yml"                  , emit: versions

    script: // see nf-core-callingcards/bin/mammals_barcode_qc.py
    """
        mammal_barcode_qc.py $bed $barcode_components \
                             $barcode_components_indicies $insert_seq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip freeze | grep pandas | sed 's/pysam==//g')
    END_VERSIONS

    """
}
