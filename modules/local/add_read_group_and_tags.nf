process ADD_RG_AND_TAGS {

    tag "$meta.id"
    label "process_high_cpu_low_mem"

    conda (params.enable_conda ? "bioconda::pysam=0.17.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam%3A0.17.0--py39h20405f9_1' :
        'quay.io/biocontainers/pysam:0.19.1--py36hea1697a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta   // fasta file
    tuple val(empty_meta), path(fai)
    val barcode_length // numeric value, eg 13
    val insertion_length



    output:
    tuple val(meta), path("*tagged.bam") , emit: bam
    path("*tagged.bam.bai")              , emit: bai
    path  "versions.yml"                 , emit: versions

    script: // see nf-core-callingcards/bin/add_read_group_to_bam.py
    """
        add_read_group_and_tags.py \\
          $bam \\
          "${meta.id}_tagged.bam" \\
          $fasta \\
          $fai \\
          $barcode_length \\
          $insertion_length

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(pip freeze | grep pysam | sed 's/pysam==//g')
    END_VERSIONS

    """
}
