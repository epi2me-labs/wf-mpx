These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore Nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository into the assets folder of Nextflow and provide a list of all parameters available for the workflow:

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
