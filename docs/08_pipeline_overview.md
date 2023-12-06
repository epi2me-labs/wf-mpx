Using community-develped tools, this workflow:
* Maps reads to a reference genome (if starting from FASTQ) (`minimap2`)
* Assesses coverage
* Discards reads not mapping to the chosen reference
* Calls variants with respect to the reference (`medaka`)
* Filters variants with <20x coverage
* Creates a draft consensus (`bcftools`)
* Creats a de-novo assembly (`flye` & `medaka`)

More information can be found in this [blog post](https://labs.epi2me.io/basic-mpox-workflow).
