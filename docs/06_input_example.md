This workflow accepts either FASTQ or BAM files as input.

The `--fastq` and `--bam` input parameters for this workflow accept a path to a single FASTQ/BAM file or a folder containing multiple FASTQ/BAM files for the sample. A sample name can be supplied with `--sample`. When BAM files are provided, these must contain alignments against the target reference.

Examples for the two possible input configurations are shown below (using FASTQ files, but the same structure can be used with BAM).

```
(i)                             (ii)
input_reads.fastq.gz        -── input_directory
                                ├── reads0.fastq.gz
                                └── reads1.fastq.gz
```