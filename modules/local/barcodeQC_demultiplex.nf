process BARCODE_QC {

    tag "$bed_meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::pysam=0.17.0 pandas=1.4.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'quay.io/biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(bed_meta), path(bed)
    tuple val(bc_meta), path(barcode_details)

    when:
    bed_meta.id == bc_meta.id

    output:
    tuple val(bed_meta), path("*_bc_fltr.bed"), emit: bed   // [val(meta), path(bed) ]
    path("*_tally.tsv" )                  , emit: tally // [ path(tally)]
    path("*_ppm.tsv")                     , emit: ppm   // [ path(ppm)]
    path  "versions.yml"                  , emit: versions

    script: // see nf-core-callingcards/bin/mammals_barcode_qc.py
    """
        barcodeQC_demultiplex.py ${bed} ${barcode_details}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip freeze | grep pandas | sed 's/pysam==//g')
    END_VERSIONS

    """
}
