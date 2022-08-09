process COUNT_HOPS {

    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.9.13 bioconda::pysam=0.17.0 conda-forge::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'quay.io/biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(bed)
    path bai
    val min_mapq

    output:
    tuple val(meta), path("*passing.bed"), path("*failing.bed"), emit: bed
    path  "versions.yml"                                       , emit: versions

    script: // see nf-core-callingcards/bin/count_hops.py
    """
        count_hops.py ${bed} ${params.require_exact_length} ${min_mapq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(pip freeze | grep pysam | sed 's/pysam==//g')
        pandas: \$(pip freeze | grep pandas | sed 's/pandas==//g')
    END_VERSIONS

    """
}
