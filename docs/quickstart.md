## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

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
