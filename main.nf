#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'

process alignReads {
    label "wfmpx"
    cpus params.align_threads
    input:
        tuple val(sample_id), val(type), path(sample_fastq), path(sample_stats)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: alignment
        tuple path("${sample_id}.bamstats"), path("${sample_id}.bam.summary"), emit: bamstats
    """
    mini_align -i ${sample_fastq} -r ${reference} -p ${sample_id} -t ${params.align_threads} -m
    stats_from_bam -o ${sample_id}.bamstats -s ${sample_id}.bam.summary -t ${params.align_threads} ${sample_id}.bam

    #Keep only mapped reads going forward
    #samtools view -b -F 4 ${sample_id}_all.bam > ${sample_id}.bam
    """
}

process coverageCalc {
    label "wfmpx"
    cpus 1
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path reference
    output:
        path "${sample_id}.depth.txt"
    """
    coverage_from_bam -s 1 -p ${sample_id} ${sample_id}.bam
    mv ${sample_id}*.depth.txt ${sample_id}.depth.txt
    """
}



process medakaVariants {
    label "wfmpx"
    cpus 2
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path reference
        path genbank
    output:
        tuple val(sample_id), val(type), path("${sample_id}.annotate.filtered.vcf")
    """
    medaka consensus ${sample_id}.bam ${sample_id}.hdf
    medaka variant --gvcf ${reference} ${sample_id}.hdf ${sample_id}.vcf --verbose
    medaka tools annotate --debug --pad 25 ${sample_id}.vcf ${reference} ${sample_id}.bam ${sample_id}.annotate.vcf

    bcftools filter -e "ALT='.'" ${sample_id}.annotate.vcf | bcftools filter -o ${sample_id}.annotate.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -

    vcf-annotator ${sample_id}.annotate.filtered.vcf ${genbank} > ${sample_id}.vcf-annotator.vcf
    """

}

process makeConsensus {
    label "wfmpx"
    cpus 1
    input:
        tuple val(sample_id), val(type), path("${sample_id}.annotate.filtered.vcf")
        path reference
        path depth
    output:
        tuple val(sample_id), val(type), path("${sample_id}.draft.consensus.fasta")
    """
    reference_name=`basename ${reference} .fasta`
    awk -v ref=\${reference_name} '{if (\$2<${params.min_coverage}) print ref"\t"\$1+1}' ${depth}  > mask.regions
    bgzip ${sample_id}.annotate.filtered.vcf
    tabix ${sample_id}.annotate.filtered.vcf.gz
    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref ${reference} -o ${sample_id}.draft.consensus.fasta ${sample_id}.annotate.filtered.vcf.gz
    """
}

process flyeAssembly {
    label "wfmpx"
    cpus params.assembly_threads
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path reference
    output:
        tuple val(sample_id), val(type), path("medaka/consensus.fasta"), path("${sample_id}_assembly_mapped.bam"), path("${sample_id}_assembly.bed")
    script:
    """
    samtools bam2fq ${sample_id}.bam > ${sample_id}_restricted.fastq
    flye --nano-raw ${sample_id}_restricted.fastq  -g 197k -t ${params.assembly_threads} --meta -o flye

    # polish assembly with medaka

    medaka_consensus -i ${sample_id}_restricted.fastq -t ${params.assembly_threads} -d flye/assembly.fasta ${params.medaka_options}

    minimap2 -ax map-ont ${reference} medaka/consensus.fasta | samtools view -bh - | samtools sort - > ${sample_id}_assembly_mapped.bam

    # make a bed file - we'll use this to plot the assembly in the report
    bedtools bamtobed -i ${sample_id}_assembly_mapped.bam > ${sample_id}_assembly.bed
    """
}

process getVersions {
    label "wfmpx"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    medaka --version | sed 's/ /,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    flye --version | sed 's/^/flye,/' >> versions.txt
    bcftools --version | head -1 | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "wfmpx"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label "wfmpx"
    input:
        path "per_read_stats/?.gz"
        path "versions/*"
        path "params.json"
        path variants
        path coverage
        path reference
        tuple val(sample_id), val(type), path(consensus), path(mapped_consensus), path(assembly_bed)
    output:
        path "wf-mpx-*.html"
    script:
        report_name = "wf-mpx-report.html"
    """

    workflow-glue report $report_name \
        --versions versions \
        per_read_stats/* \
        --params params.json \
        --variants ${variants} \
        --coverage ${coverage} \
        --reference ${reference} \
        --assembly_bed ${assembly_bed}
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfmpx"
    publishDir (
        "${params.out_dir}",
        mode: 'copy',
        saveAs: { fname_as ?: fname }
    )
    input:
        tuple path(fname), val(fname_as)
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        reference
        genbank
    main:
        samples_for_processing = reads.map {it -> [it[0].alias, it[0].type, it[1], it[2]]}
        alignment = alignReads(samples_for_processing, reference)
        coverage = coverageCalc(alignment.alignment, reference)
        variants = medakaVariants(alignment.alignment, reference, genbank)
        draft = makeConsensus(variants, reference, coverage)
        software_versions = getVersions()
        workflow_params = getParams()

        if (params.assembly == true) {
          assembly = flyeAssembly(alignment.alignment, reference)
        }


        report = makeReport(
            reads.map{ it[2].resolve("per-read-stats.tsv.gz") }.toList(),
            software_versions.collect(),
            workflow_params,
            variants.map{it[2]}.collect(),
            coverage,
            reference,
            assembly
        )
    emit:
        results = report.concat(
            alignment.alignment.map{it[2]}.collect(),
            alignment.alignment.map{it[3]}.collect(),
            variants.map{it[2]}.collect(),
            variants.map{it[3]}.collect(),
            draft.map{it[2]}.collect(),
            coverage,
            assembly.map{it[2]}.collect(),
            software_versions
        )
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.reference == null){
        params.remove('reference')
        params._reference = projectDir.resolve("./data/references/MT903344.1.fasta").toString()
        params._genbank = projectDir.resolve("./data/references/MT903344.1.gb").toString()
    } else {
        path = projectDir.resolve("./data/references/"+params.reference+".fasta").toString()
        params._reference = file(path, type: "file", checkIfExists:true).toString()
        path_genbank = projectDir.resolve("./data/references/"+params.reference+".gb").toString()
        params._genbank = file(path_genbank, type: "file", checkIfExists:true).toString()
        params.remove('reference')
    }

    samples = fastq_ingress([
       "input":params.fastq,
       "sample":params.sample,
       "sample_sheet":null,
       "stats":true ])

    pipeline(samples, params._reference, params._genbank)

    pipeline.out.results.map {
        // map outputs to [it, null] to write outputs to top level out_dir
        it -> [it, null]
    }.mix(
        // mix in per-read-stats
        samples.map {it -> [it[2].resolve("per-read-stats.tsv.gz"), "${it[0].alias}.per-read-stats.tsv.gz"]}
    ) 
    | output

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
