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

## DEVELOPMENT INSTALLATION

You will need the following two pieces of software to run this pipeline:

1. [Nextflow](https://www.nextflow.io/)
2. One of: [Singularity](https://sylabs.io/singularity/), [Docker](https://www.docker.com/) or [conda](https://docs.conda.io/en/latest/) (singularity or docker are far preferred)

AND, you will need to download [git lfs](https://git-lfs.github.com/) in order to
fully clone this repository. This is not a good way of doing the test data, but
until I curtail the chr1 fasta file, it is how I am doing it.
If you are on a linux system that you own, do this:

```
$ sudo apt update
$ sudo apt install git-lfs
$ git lfs install
```
On the HTCF cluster, you can do this through `spack`.

If you're using mac, then you can probably substitute `apt` for `brew`. If you're
using Windows, you should think about the life choices that have brought you to
this point, and then partition your hard drive and install ubuntu. Or, google
about and then let me know what works, and I'll update the readme.

Alternatively, email me, and I'll send you a a zipped version of the repo.

Note: it isn't necessary to make the cc_tester directory -- call it whatever you
like, or do this differently. Just make sure the path in the run script below
points to the main.nf file in the callingcards repo.

Assuming you have git lfs, then clone the repo.

```
$ mkdir cc_tester
# you could download/install the nextflow executable in this directory
$ cd cc_tester
$ git clone https://github.com/cmatKhan/callingcards.git
```

Next, copy and paste the script below into a file called, for example, `run_nf.sh`

```
#!/bin/bash

mkdir tmp

# CHOOSE ONE OF SINGULARITY, DOCKER OR CONDA. I haven't tested docker. It probably works
nextflow run callingcards/main.nf  -profile test_human,<singularity/docker/conda> -resume

```
Note that this assumes that nextflow is either in your `$PATH`, or the executable is
in the same directory from which you are launching this script.

Make the script executable (on a linux system, `chmod +x run_nf.sh`) and then launch:

```
./run_nf.sh
```

You could just execute the line in the runscript, but if there are errors, I'll ask to see the
submission command -- saving the run command is useful for debugging.

This will run the pipeline locally, meaning the compute resources come from
the machine from which you launch. Be careful -- don't do this from a login node
on the cluster, for example. This will work on your local computer, or an interactive
node, with at least 16GB ram.

If you want to test this via scheduler submission, then here are examples of
SLURM and and SGE submission scripts:

__SLURM__ (htcf)
```
#!/usr/bin/env bash

eval $(spack load --sh openjdk)
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow)

mkdir tmp

nextflow run callingcards/main.nf  -profile test_slurm,singularity -resume
```


__SGE__ (dsg -- the genetics cluster. Note that I'm using conda, not singularity, since singularity isn't available on their cluster)
```
#!/bin/bash

#$ -N nf_tester
#$ -cwd
#$ -o /dsgmnt/llfs_external/$USER/logs
#$ -e /dsgmnt/llfs_external/$USER/logs

root=/dsgmnt/llfs_external/cmateusiak

# activate environment
# note that in this case, I installed nextflow using conda like so:
# conda create -p /path/to/envs/dir/nextflow nextflow
source activate $root/conda_envs/nextflow

mkdir tmp

nextflow run callingcards/main.nf  -profile test_sge,conda -resume

```

### yeast data

A test profile for yeast hasn't yet been created, but this will handle yeast (paired or single end data, with barcodes either only on the first read, or on both reads) data through the pile up, also. Here is an example of what a params file, and the input sheet, would look like:

params.json
```
{
  "input":"input_samplesheet.csv",
  "aligner":"bwamem2",
  "fasta":"\/home\/oguzkhan\/ref\/S288C_R64\/GCF_000146045.2_R64_genomic.fa",
  "r1_bc_pattern":"NNNNNXXXXXXXXXXXXXXXXX",
  "r2_bc_pattern":"NNNNNNNNXXXX",
  "barcode_length": 13,
  "samtools_mpileup_adjust_mq": 50,
  "min_mapq": 10,
  "promoter_bed": "\/home\/oguzkhan\/Desktop\/tmp\/cc_tester\/rmlab_notOrf.bed",
  "background_data": "\/home\/oguzkhan\/Desktop\/tmp\/cc_tester\/old_testin\/NOTF_Minus_Adh1_2015_17_combined_chase_edit.csv",
  "pileup_stranded": "FALSE"
}

```

input_sheet.csv
```
sample,fastq_1,fastq_2,barcodes
run_5399,run_5399/PhiX_S1_R1_001.fastq.gz,run_5399/PhiX_S1_R2_001.fastq.gz,/home/oguzkhan/Desktop/tmp/cc_tester/old_testing/run_5399_demult_barcodes.tsv

```
And this is what demult barcode looks like for this run:

```
MIG2	TCAGTCCCGTTGG
CAT8	GCCTGGGCGGCAG
GLN3	ATTTGGGGGGGGT
ARO80	TTGGTGGGGGTAG
CBF1	CTCGGTCGTCAGT
```
If you want to test the pipeline on yeast data, it is likely best just to let me know and let me help get it started (until I get a test profile made, at which point, it will be like above, where you can just use the profile 'test_yeast' rather than 'test_human')

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
