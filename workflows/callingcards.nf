/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCallingcards.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input,
                           params.multiqc_config,
                           params.fasta,
                           params.fasta_index,
                           params.bwamem2_index ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK            } from '../subworkflows/local/1_input_check'
include { SAMTOOLS_INDEX_FASTA  } from '../subworkflows/nf-core/2_samtools_index_fasta'
include { UMITOOLS_FASTQC        } from '../subworkflows/nf-core/3_umitools_fastqc'
include { ALIGN                  } from '../subworkflows/local/4_align'
include { PROCESS_ALIGNMENTS     } from '../subworkflows/local/5_process_alignments'
include { QUANTIFY_HOPS          } from '../subworkflows/local/6_quantify_hops'
include { PROCESS_QUANTIFICATION } from '../subworkflows/local/7_process_quantification'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
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

    //
    // SUBWORKFLOW_1: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW_2: use samtools to create a .fai index for the genome
    // input:
    // output:
    //
    // if the user does not provide an genome index, index it
    if (!params.fasta_index){
        SAMTOOLS_INDEX_FASTA ( params.fasta )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FASTA.out.versions.first())
        ch_fasta_index = SAMTOOLS_INDEX_FASTA.out.fai
    } else {
        ch_fasta_index = Channel.of(["",params.fasta_index])
    }

    //
    // SUBWORKFLOW_3: run sequencer level QC, extract barcodes and trim
    //
    UMITOOLS_FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(UMITOOLS_FASTQC.out.versions.first())

    //
    // SUBWORKFLOW_4: align reads
    // input:
    // output:
    //
    ALIGN (
        UMITOOLS_FASTQC.out.reads,
        params.fasta
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions.first())

    //
    // SUBWORKFLOW_5: sort, add barcodes as read group, add tags, index and
    //              extract basic alignment stats
    //
    PROCESS_ALIGNMENTS (
        ALIGN.out.bam,
        params.fasta,
        ch_fasta_index
    )
    ch_samtools_stats    = PROCESS_ALIGNMENTS.out.stats
    ch_samtools_flagstat = PROCESS_ALIGNMENTS.out.flagstat
    ch_samtools_idxstats = PROCESS_ALIGNMENTS.out.idxstats
    ch_versions          = ch_versions.mix(PROCESS_ALIGNMENTS.out.versions.first())

    //
    // SUBWORKFLOW_6: turn alignments into ccf (modified bed format) which
    //              may be used to quantify hops per TF per promoter region
    // QUANTIFY_HOPS (
    //     PROCESS_ALIGNMENTS.out.bed
    //     barcode_length,
    //     promoter_bed,
    //     background_data,
    //     pileup_stranded
    // )

    //
    // SUBWORKFLOW_7: calculate statistics and other metrics relating to the hops
    //
    // PROCESS_QUANTIFICATION (

    // )

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
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_FASTQC.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))
    // add BAMQC output here, also -- see RNAseq pipeline


    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
