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


