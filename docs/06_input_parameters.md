### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bam | string | BAM file to use in the analysis. | Path to a single BAM file, or a single directory containing one or more BAM files. Multiple samples are not currently supported. |  |
| fastq | string | FASTQ file to use in the analysis. | Path to a single FASTQ file, or a single directory containing one or more FASTQ files. Multiple samples or barcodes are not currently supported. |  |
| reference | string | The reference genome to use for mapping. | This is used if inputting a FASTQ file. We provide five MPXV reference sequences for alignment, representing clades I and II. More information can be found in our [blog post introducing this workflow](https://labs.epi2me.io/basic-mpox-workflow) | NC_003310.1 |


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


