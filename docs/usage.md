# nf-core/callingcards: Usage

<!-- ## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/callingcards/usage](https://nf-co.re/callingcards/usage) -->

<!-- > _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._ -->

## Introduction

Calling Cards experiments may be performed in both yeast and mammalian cells. The processing steps diverge, and this
divergence is controlled by setting certain parameters. Suggested default parameters for yeast and mammalian
processing runs are provided through the profiles [yeast](../conf/default_yeast.config) and [mammal](../conf/default_mammal.config). These may be used by simply including them with the `-profile` flag, for instance:

```
$ nextflow run callingcards/main.nf \
    -profile default_yeast,singularity \
    -c local.config \
    --input /path/to/samplesheet.csv
    --fasta /path/to/genome.fasta
    --output results_20220811
```

`local.config` is a [configuration](https://nf-co.re/usage/configuration) file which is the best place to put configuration settings such as what [executor](https://www.nextflow.io/docs/latest/executor.html)
you wish to use. If your institution already has a configuration profile on nf-core, then you should use that profile in the `-profile` flag instead.nf-

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```
Or, if you are using a file to save the run parameters (recommended), rather than submitting them on the command line,
then the file, with only the input set, would look like so:

```json

{
  "input":"input_samplesheet.csv"
}

```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2,barcode_details
mouse_AY60-6_50,mouse/test_data/AY60-6_50k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY60-6_100,mouse/test_data/AY60-6_100k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_50_lowQuality,mouse/test_data/AY09-1_50k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_100_lowQuality,mouse/test_data/AY09-1_100k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json

```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`          | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`         | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`         | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".
| `barcode_details` | Full path to the barcode details json file for a given sample. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Barcode Details

The barcode details json stores data which allows the pipeline to relate sequence barcodes in the calling cards reads to a given transcription factor.

#### [Yeast](../assets/yeast/barcode_details.json)


```json

{
    "indicies": {
        "r1_primer_bc":     [0, 5],
        "transposon_seq":   [5, 22],
        "r2_primer_bc":     [22, 30],
        "restriction_site": [30, 34]
    },
    "components": {
        "r1_primer_bc":     ["TCAGT", "GCCTG", "ATTTG", "TTGGT", "CTCGG"],
        "transposon_seq":   ["AATTCACTACGTCAACA"],
        "r2_primer_bc":     ["CCCGTTGG", "GGCGGCAG", "GGGGGGGT", "GGGGGTAG", "TCGTCAGT"],
        "restriction_site": ["TCGA","GCGC","CCGG"]
    },
    "tf_map": {
        "r1_primer_bc":    ["TCAGT", "GCCTG", "ATTTG", "TTGGT", "CTCGG"],
        "r2_primer_bc":    ["CCCGTTGG", "GGCGGCAG", "GGGGGGGT", "GGGGGTAG", "TCGTCAGT"],
        "TF":              ["MIG2", "CAT8", "GLN3", "ARO80", "CBF1"]
    }

}

```

#### [Mammals](../assets/human/barcode_details.json)

```json

{
    "indicies": {
        "om_pb": [0, 3],
        "pb_lrt1": [3, 28],
        "srt": [28, 32],
        "pb_lrt2": [32, 38]
    },
    "components": {
        "om_pb":   ["TAG"],
        "pb_lrt1": ["CGTCAATTTTACGCAGACTATCTTT"],
        "srt":     ["CTAG", "CAAC", "CTGA", "GCAT", "GTAC", "CACA", "TGAC", "GTCA",
                    "CGAT", "CTCT", "GAAG", "TCGA", "CATG", "GTTG", "CTTC", "GCTA",
                    "GAGA", "GTGT", "CGTA", "TGGT", "GGAA", "ACAC", "TCAG", "TTGG",
                    "CAGT", "TTTT"],
        "pb_lrt2": ["GGTTAA"]
    },
    "insert_seq": ["TTAA"]

}

```

## Running the pipeline

The typical command for running the pipeline is as follows:

```
$ nextflow run callingcards/main.nf \
    -profile default_mammals,singularity \
    -c local.config \
    -params-file params.json
```

This will launch the pipeline with the `default_mammal` pipeline settings. It will use the `singularity` profile, and
set user specific system configuration (eg executor settings) in a file called local.config. Parameters such as `input`,
the path to the samplesheet, `output`, and other run parameters will be set in the params.json file.
See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Running the pipeline on HTCF

You will need the following system dependencies. These can be installed with
`spack`, and therefore only need to be done one time by one person in the lab.

```

$ interactive

$ spack install git

$ spack install openjdk

$ spack install nextflow

$ spack install singularityce

```

Next, navigate into your scratch space, make a directory, and pull the pipeline
repo if you haven't already:

```
$ cd /scratch/<lab>/<your scratch folder>

# it doesn't matter what this is called, just that you know what it is
$ mkdir nf_cc

$ cd nf_cc

# if someone hasn't installed git with spack in your lab yet, then launch an
# interactive session and install git first
$ eval $(spack load --sh git)

$ git clone https://github.com/cmatKhan/callingcards
```

Copy and paste the script below into a file called, for example, `run_nf.sh`

```bash
#!/usr/bin/env bash

#SBATCH --mem-per-cpu=10G
#SBATCH -J cc_nf_test.out
#SBATCH -o cc_nf_test.out

# load system dependencies -- on HTCF, we use spack
eval $(spack load --sh openjdk)
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow)

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir singularity

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

# if the repo is in this directory, then this relative path callingcards/main.nf
# will work. Otherwise, replace callingcards/main.nf with the path, relative or
# absolute, to the correct main.nf
nextflow run callingcards/main.nf -profile test_slurm,singularity -resume
```

Launch this with the command

```bash
$ sbatch run_nf.sh

```
At this point, you can check progress with `squeue -u $USER`. You can also
track progress by looking at the nextflow process log, which in this case
will be called `cc_nf_test.out`:

```bash
$ tail -150 cc_nf_test.out
```
Note that it sometimes takes HTCF a long time to start scheduling processes. If
there is no error message in `cc_nf_test.out`, then just keep waiting.

### Running the pipeline on your laptop

See the [README]("../README.md") Quick Start section for instructions on installing
nextflow and one of docker, singularity or conda (conda should be an absolute last choice).
Then follow the instructions in the HTCF tutorial above to create a directory
in which you will clone the callingcards pipeline repository, and launch the pipeline.
Next, copy and paste the script below into a script named something like `run_nf.sh`

```bash

#!/usr/bin/bash

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)
mkdir singularity

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

nextflow run nf-core-callingcards/main.nf  -profile test_yeast,singularity -resume

```
You'll run this with the command `.run_nf.sh` and the process logger will
print to `stdout`.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
$ cd /repo/path/callingcards

$ git pull
```

#### note this won't work unless this goes onto nf-core
```console
nextflow pull nf-core/callingcards
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/callingcards releases page](https://github.com/nf-core/callingcards/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
