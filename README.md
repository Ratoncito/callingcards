# ![nf-core/callingcards](docs/images/nf-core-callingcards_logo_light.png#gh-light-mode-only) ![nf-core/callingcards](docs/images/nf-core-callingcards_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/callingcards/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/callingcards/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/callingcards/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/callingcards/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/callingcards/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23callingcards-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/callingcards)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/callingcards** is a bioinformatics best-practice analysis pipeline for A bioinformatics analysis pipeline for processing Transposon Calling Cards sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

The testing data is currently provided in the `assets` directory of the pipeline code repository. Instructions on running the tests area available in Usage.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/callingcards/results).

## Pipeline summary

1. Prepare Reads
    1. Extract barcodes ([UMItools](https://github.com/CGATOxford/UMI-tools))
    1. Trim, and reduce to only R1 depending on user input ([Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic))
1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
1. Prepare the Genome
    1. Samtools faidx and aligner indicies
1. Alignment
    1. [bwamem2](https://github.com/bwa-mem2/bwa-mem2)
1. Process Alignments
    1. Extract alignment QC metrics ([samtools](https://www.htslib.org/))
    1. Add Read Group and tags to alignment files ([custom script](https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py))
    1. Quantify transposon hops ([custom script](https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py))
1. Process Transposon Quantification
    1. Demultiplex Transcription Factors ([custom script](https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py))
    1. Calling Cards specific quality assessment ([custom script](https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py))
    1. Call peaks and calculate significance ([custom script](https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py))
1. Present QC for raw read and alignment metrics ([`MultiQC`](http://multiqc.info/))

## Quick Start

Note that more detailed instructions are available in [usage](https://nf-co.re/callingcards/usage), including examples of launching the Calling Cards pipeline in specific environments. There are currently detailed examples for launching the pipeline on a laptop, and a cluster which runs [SLURM](https://slurm.schedmd.com/overview.html) and manages dependencies with [spack](https://spack.readthedocs.io/en/latest/)).

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

1. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

1. Git clone the repo and test it on a minimal data set designed to confirm that the pipeline will complete

   ```bash
   $ mkdir calling_card_output

   $ cd calling_card_output

   $ git clone https://github.com/cmatKhan/callingcards/blob/main/bin/barcodeQC_demultiplex.py

   # you can replace singularity below with any one of the dependency managers
   # There are also other test profiles -- see usage for more detail
   $ nextflow run callingcards/main.nf -profile test_yeast,singularity
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`test_yeast` and `singularity` in the example command above). You can chain multiple config profiles in a comma-separated string, as demonstrated.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

1. Start running your own analysis!

   ```bash
   nextflow run callingcards/main.nf -params-file params.json -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

The `params.json` file is described in [usage](https://nf-co.re/callingcards/usage)

## Documentation

The nf-core/callingcards pipeline comes with documentation about the pipeline [usage](https://nf-co.re/callingcards/usage), [parameters](https://nf-co.re/callingcards/parameters) and [output](https://nf-co.re/callingcards/output).

## Credits

nf-core/callingcards was originally written by Chase Mateusiak.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#callingcards` channel](https://nfcore.slack.com/channels/callingcards) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/callingcards for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


## DEVELOPMENT INSTALLATION

### Pipeline Dependencies

You will need the following two pieces of software to run this pipeline:

1. [Nextflow](https://www.nextflow.io/)
2. One of: [Singularity](https://sylabs.io/singularity/),
[Docker](https://www.docker.com/) or
[conda](https://docs.conda.io/en/latest/) (singularity or docker are far preferred)

### SETUP

I suggest creating a directory from which you will launch the pipeline. The
name of this directory doesn't matter.

```bash
# make a new directory in your current working directory
$ mkdir cc_tester
# change locations into the new directory
$ cd cc_tester
```
Next, clone this repository

```bash
$ git clone https://github.com/cmatKhan/callingcards.git
```

Below are examples scripts used to launch the pipeline on a local computer,
or on HTCF (using slurm). Copy and paste the script into a file called,
for instance, `run_nf.sh`. Note that the name of this file does not matter.

### Running the pipeline in different environments

__LOCAL__ (your laptop)

```
#!/bin/bash


tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir singularity

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp


# CHOOSE ONE OF SINGULARITY, DOCKER OR CONDA. Conda itself is buggy -- use only as a very last resort.
# CHOOSE EITHER test_human or test_yeast
nextflow run callingcards/main.nf  -profile <test_human, test_yeast>,<singularity/docker/conda> -resume

```
Make the script executable (on a linux system, `chmod +x run_nf.sh`) and then launch:

Note that this assumes that nextflow is either in your `$PATH`, or the executable is
in the same directory from which you are launching this script.

To launch the script, do the following:

```
$ ./run_nf.sh
```

__SLURM__ (htcf)
```
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

nextflow run callingcards/main.nf -profile test_slurm,singularity -resume
```


# IGNORE EVERYTHING BELOW
## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/callingcards** is a bioinformatics analysis pipeline for processing Transposon Calling Cards sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/callingcards/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Check the sample sheet
2. Append barcodes to the fastq ID line with UMI Tools ([`UMI Tools`](https://umi-tools.readthedocs.io/en/latest/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align ([`bwamem2 mem`](https://github.com/bwa-mem2/bwa-mem2https://umi-tools.readthedocs.io/en/latest/))
4. Add the barcodes as read groups to the bam, index and sort. Implemented in ([a custom script](https://github.com/BrentLab/callingcards/blob/main/bin/add_read_group.py) using [pysam](https://pysam.readthedocs.io/en/latest/api.html).
5. Bam QC ([`Qualimap`](http://qualimap.conesalab.org/))
6.
5. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/callingcards -profile test,YOURPROFILE
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```console
    nextflow run nf-core/callingcards -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

## Documentation

The nf-core/callingcards pipeline comes with documentation about the pipeline [usage](https://nf-co.re/callingcards/usage), [parameters](https://nf-co.re/callingcards/parameters) and [output](https://nf-co.re/callingcards/output).

## Credits

nf-core/callingcards was originally written by Chase Mateusiak, Woo Jung.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#callingcards` channel](https://nfcore.slack.com/channels/callingcards) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/callingcards for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
