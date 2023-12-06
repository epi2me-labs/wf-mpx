# Mpox metagenomics assembly workflow

A basic workflow for analysing ONT Mpox data to create draft consensus or de-novo assemblies.



## Introduction

This workflow provides a simple way to analyse mpox sequencing data; taking
raw Oxford Nanopore Technologies reads and creating a draft consensus and assembly.

We expect data to be generated from a mixed population sample hence the use of the term "metagenomics".

> No trimming of sequences is carried out so be vigilant when using targeted data.




## Compute requirements

Recommended requirements:

+ CPUs = 4
+ Memory = 16GB

Minimum requirements:

+ CPUs = 4
+ Memory = 10GB

Approximate run time: 45 minutes

ARM processor support: True




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow:

```
nextflow run epi2me-labs/wf-mpx -help 
```

A demo dataset is provided for testing of the workflow. It can be downloaded using: 
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-mpx/wf-mpx-demo.tar.gz
tar -xzvf wf-mpx-demo.tar.gz
```
The demo data can then be run with this command:
```
nextflow run epi2me-labs/wf-mpx --fastq wf-mpx-demo/fastq/barcode01/barcode01.fastq.gz
```

Either FASTQ or BAM files can be used as input.




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Inputs

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bam | string | BAM file to use in the analysis. | Path to a single BAM file, or a directory containing a single BAM file. Multiple samples are not currently supported. |  |
| fastq | string | FASTQ file to use in the analysis. | Path to a single FASTQ file, or a directory containing a single FASTQ file. Multiple samples are not currently supported. |  |
| reference | string | The reference genome to use for mapping. | This is used if inputting a FASTQ file. We provide four popular mpox reference sequences with which to map your reads to. More information can be found here: https://labs.epi2me.io/basic-mpox-workflow | MT903344.1 |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| assembly | boolean | Perform assembly with flye. | Take the reads mapped to the mpx reference and perform a de novo assembly using flye. This might be useful to uncover any structural variation not present in your chosen reference genome. Turning off will reduce run times and computational intensity at the expense of not being able to identify structural variants. | True |
| min_coverage | number | Coverage threshold for masking in the consensus step. | Regions with less than the coverage entered here are masked (nucleotide sequences will be replaced with N) when generating the consensus. It might be useful to decrease this value if you have very low coverage, but this will affect the quality of your consensus sequence. | 20 |
| medaka_options | string | Pass through options to `medaka_consensus` during assembly. | Only applicable when flye assembly is used (default). You can use this field to give `medaka_consensus` command line options for the assembly process of the workflow. For example, you could change the model by entering `"-m r941_min_high_g303"`. The full string of options should be quoted as in the example. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Outputs files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-mpx-report.html | The report for the workflow | aggregated |
| Assembly FASTQ | ./{{ alias }}.final.fastq | Sequence and quality score for final assembly. | per-sample |
| Consensus assembly FASTA | ./consensus.fasta | De-novo consensus assembly sequence from flye and polished by medaka. | per-sample |
| Draft consensus FASTA | ./{{ alias }}.draft.consensus.fasta | Draft consensus sequence from bcftools. | per-sample |
| Read Stats | ./{{ alias }}.per-read-stats.tsv.gz | A simple text file providing a summary of sequencing reads. | per-sample |
| Read alignment | ./{{ alias }}.bam | Read alignments in BAM format. | per-sample |
| Variants file | ./{{ alias }}.annotate.filtered.vcf | Called variants in VCF format. | per-sample |




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

## FAQs
If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-mpx/issues) page or start a discussion on the [community](https://nanoporetech.com/support).




## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




