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
