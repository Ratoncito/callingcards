/*
================================================================================
    VALIDATE INPUTS
================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCallingcards.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input,
                           params.multiqc_config,
                        //    params.fasta,
                           params.fasta_index,
                           params.bwamem2_index,
                           params.promoter_bed,
                           params.background_qbed,
                           params.chr_map ]

for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
    }
}

// Check mandatory parameters
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}

/*
================================================================================
    CONFIG FILES
================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml",
                            checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ?
                            Channel.fromPath(params.multiqc_config) :
                            Channel.empty()

/*
================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from "${projectDir}/subworkflows/local/1_input_check"
include { PREPARE_GENOME     } from "${projectDir}/subworkflows/local/2_prepare_genome"
include { PREPARE_READS      } from "${projectDir}/subworkflows/local/3_prepare_reads"
include { ALIGN              } from "${projectDir}/subworkflows/local/4_align"
include { PROCESS_ALIGNMENTS } from "${projectDir}/subworkflows/local/5_process_alignments"
include { PROCESS_HOPS       } from "${projectDir}/subworkflows/local/6_process_hops"
include { STATISTICS         } from "${projectDir}/subworkflows/local/7_statistics"

/*
================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
================================================================================
    RUN MAIN WORKFLOW
================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CALLINGCARDS {

    // instantiate channels
    ch_versions          = Channel.empty()
    ch_fasta_index       = Channel.empty()
    ch_bam_index         = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flatstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()

    fasta = params.fasta?
            Channel.fromPath(params.fasta).collect():
            Channel.empty()

    ch_promoter_bed      = Channel.of(params.promoter_bed).collect()
    ch_background_qbed   = Channel.of(params.background_qbed).collect()
    ch_chr_map           = Channel.of(params.chr_map).collect()

    //
    // SUBWORKFLOW_1: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW_2: Index the genome in various ways
    // input:
    // output:
    //
    // if the user does not provide an genome index, index it
    PREPARE_GENOME(
        fasta
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW_3: run sequencer level QC, extract barcodes and trim
    //
    PREPARE_READS (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(PREPARE_READS.out.versions)

    //
    // SUBWORKFLOW_4: align reads
    // input:
    // output:
    //
    ALIGN (
        PREPARE_READS.out.reads,
        PREPARE_GENOME.out.bwamem2_index
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    //
    // SUBWORKFLOW_5: sort, add barcodes as read group, add tags, index and
    //              extract basic alignment stats
    //
    PROCESS_ALIGNMENTS (
        ALIGN.out.bam,
        fasta,
        PREPARE_GENOME.out.fai
    )
    ch_versions          = ch_versions.mix(PROCESS_ALIGNMENTS.out.versions)

    // join the barcode details file back to the corresponding samples
    PROCESS_ALIGNMENTS.out.bed
        .map{ meta, passing, failing -> [meta, [passing, failing]]}
        .mix(INPUT_CHECK.out.barcode_details)
        .groupTuple()
        .multiMap{meta, file_list ->
            a: [ add_status_to_meta(meta,file_list[1][0]), file_list[0], file_list[1][0] ]
            b: [ add_status_to_meta(meta,file_list[1][1]), file_list[0], file_list[1][1] ] }
        .set{ ch_passing_failing_bed }

    ch_passing_failing_bed.a
        .mix(ch_passing_failing_bed.b)
        .set{ process_hops_input }

    // SUBWORKFLOW_6a: turn alignments into qbed (modified bed format) which
    //              may be used to quantify hops per TF per promoter region
    PROCESS_HOPS (
        process_hops_input
    )
    ch_versions = ch_versions.mix(PROCESS_HOPS.out.versions)

    STATISTICS (
        PROCESS_HOPS.out.qbed,
        ch_promoter_bed,
        ch_background_qbed,
        ch_chr_map,
        params.standard_chr_format,
        params.sqlite_db_out,
        params.poisson_pseudocount
    )
    ch_versions = ch_versions.mix(PROCESS_HOPS.out.versions)

    //
    // collect software versions into file
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCallingcards.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.umi_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.trimmomatic_log.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.picard_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.preseq_lcextrap.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.preseq_ccurve.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.idxstats.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
================================================================================
    COMPLETION EMAIL AND SUMMARY
================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
================================================================================
    THE END
================================================================================
*/

/*
================================================================================
    FUNCTIONS
================================================================================
*/

def add_status_to_meta(Map meta, filename) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    def (_,status) = (filename =~ /(passing|failing)/ )[0]

    new_meta["status"] = status

    return new_meta
}
