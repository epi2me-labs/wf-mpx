# MPXV metagenomic assembly workflow

A workflow for analysing ONT monkeypox virus (MPXV) sequences to create draft consensus or de novo assemblies.



## Introduction

This workflow provides a simple way to analyse mpox sequencing data; taking
raw Oxford Nanopore Technologies reads and creating a draft consensus and assembly.

We expect data to be generated from a mixed population sample hence the use of the term "metagenomics".

> No trimming of sequences is carried out so be vigilant when using targeted data.




## Compute requirements

Recommended requirements:

+ CPUs = 4
+ Memory = 32GB

Minimum requirements:

+ CPUs = 4
+ Memory = 16GB

Approximate run time: 45 minutes

ARM processor support: True




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-mpx --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-mpx
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-mpx/wf-mpx-demo.tar.gz
tar -xzvf wf-mpx-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-mpx \
	--fastq 'wf-mpx-demo/fastq/barcode01' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

This workflow accepts either FASTQ or BAM files as input.

The `--fastq` and `--bam` input parameters for this workflow accept a path to a single FASTQ/BAM file or a folder containing multiple FASTQ/BAM files for the sample. A sample name can be supplied with `--sample`. When BAM files are provided, these must contain alignments against the target reference.

Examples for the two possible input configurations are shown below (using FASTQ files, but the same structure can be used with BAM).

```
(i)                             (ii)
input_reads.fastq.gz        -── input_directory
                                ├── reads0.fastq.gz
                                └── reads1.fastq.gz
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bam | string | BAM file to use in the analysis. | Path to a single BAM file, or a single directory containing one or more BAM files. Multiple samples are not currently supported. |  |
| fastq | string | FASTQ file to use in the analysis. | Path to a single FASTQ file, or a single directory containing one or more FASTQ files. Multiple samples or barcodes are not currently supported. |  |
| reference | string | The reference genome to use for mapping. | This is used if inputting a FASTQ file. We provide five MPXV reference sequences for alignment, representing clades I and II. More information can be found in our [blog post introducing this workflow](https://labs.epi2me.io/basic-mpox-workflow) | NC_003310.1 |
| reference_fasta | string | A reference genome FASTA file to use for mapping. | This is used if inputting a FASTQ file and when wanting to use a dfferent reference than the options provided with the workflow through the '--reference' parameter. |  |
| reference_gb | string | The reference genome genbank (.gb) file to use for mapping. | This is used if inputting a FASTQ file and when wanting to use a dfferent reference than the options provided with the workflow through the '--reference' parameter. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample | string | A sample name for the outputs. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |
| assembly | boolean | Perform assembly with flye. | Take the reads mapped to the mpx reference and perform a de novo assembly using Flye. This might be useful to uncover any structural variation not present in your chosen reference genome. Turning off will reduce run times and computational intensity at the expense of not being able to identify structural variants. | True |
| min_coverage | number | Coverage threshold for masking in the consensus step. | Regions with less than the coverage entered here are masked (nucleotide sequences will be replaced with N) when generating the consensus. It might be useful to decrease this value if you have very low coverage, but this will affect the quality of your consensus sequence. | 20 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-mpx-report.html | The report for the workflow | aggregated |
| De novo consensus assembly FASTA | ./denovo.consensus.fasta | De novo consensus assembly sequence from Flye and polished by Medaka. | per-sample |
| Reference-based consensus assembly FASTA | ./{{ alias }}.ref.consensus.fasta | Reference-based consensus sequence from Bcftools. | per-sample |
| Read stats | ./{{ alias }}.per-read-stats.tsv.gz | A simple text file providing a summary of sequencing reads. | per-sample |
| Read alignment | ./{{ alias }}.bam | Read alignments in BAM format. | per-sample |
| Alignment index file | ./{{ alias }}.bam.bai | Index file of BAM file. | per-sample |
| Variants file | ./{{ alias }}.annotate.filtered.vcf | Called variants in VCF format. | per-sample |
| Depth file | ./{{ alias }}.annotate.filtered.vcf | Per-base depth: overall, forward and reverse. | per-sample |




## Pipeline overview

Using community-develped tools, this workflow:
* Maps reads to a reference genome (if starting from FASTQ) (`minimap2`)
* Assesses coverage
* Discards reads not mapping to the chosen reference
* Calls variants with respect to the reference (`medaka`)
* Filters variants with <20x coverage
* Creates a draft consensus (`bcftools`)
* Creats a de-novo assembly (`flye` & `medaka`)

More information can be found in this [blog post](https://labs.epi2me.io/basic-mpox-workflow).




## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-mpx/issues) page or start a discussion on the [community](https://nanoporetech.com/support).




## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




