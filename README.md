# wf-mpx | Mpox Virus Metagenomic Assembly

This repository contains a [nextflow](https://www.nextflow.io/) workflow for the
assembly of reads originating from the mpox virus obtained through
Oxford Nanopore metagenomic sequencing.
## Introduction

This workflow provides a simple way to analyse mpox sequencing data; taking
raw Oxford Nanopore Technologies reads and creating a draft consensus and assembly.

We expect data to be generated from a mixed population sample hence the use of the term "metagenomics".

> No trimming of sequences is carried out so be vigilant when using targeted data.

Using community-develped tools, this workflow:
* Maps reads to a reference genome (`minimap2`)
* Assesses coverage
* Discards reads not mapping to the chosen reference
* Calls variants with resepect to the reference (`medaka`)
* Filters variants with <20x coverage
* Creates a draft consensus (`bcftools`)
* Creats a de-novo assembly (`flye` & `medaka`)

More information can be found in this [blog post](https://labs.epi2me.io/basic-mpox-workflow).
## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-mpx --help
```

to see the options for the workflow.

**Workflow example**

To run the workflow with test data

```
git clone https://github.com/epi2me-labs/wf-mpx
nextflow run epi2me-labs/wf-mpx --fastq wf-mpx/test_data/fastq/barcode01
```

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow.
* a draft consensus sequence obtained
* a medaka polished assembly
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
