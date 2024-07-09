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
include { xam_ingress } from './lib/ingress'  

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process alignReads {
    label "wfmpx"
    cpus 4
    memory '2GB'
    input:
        // `reads` can be FASTQ, BAM, or uBAM
        tuple val(meta), path("input"), path("index")
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: alignment
        tuple path("${sample_id}.bamstats"), path("${sample_id}.bam.summary"), emit: bamstats
    script:
        sample_id = meta.alias
        type = meta.type
    """
    # if FASTQ or uBAM, align
    if [[ ${params.fastq || meta.is_unaligned} = true ]]; then
        ${params.fastq ? "cat" : "samtools fastq -T '*'" } input \\
        | minimap2 -ax map-ont --secondary=no -L --MD -t 3 "$reference" - \\
            --cap-kalloc 100m --cap-sw-mem 50m \\
        | samtools sort --write-index -o ${sample_id}.bam##idx##${sample_id}.bam.bai
    else
        # we got BAM; make sure it was aligned against the provided reference
        samtools faidx "$reference"
        ref_name=\$(cut -f1 ${reference}.fai)
        bam_ref_name=\$(
            samtools view -H input | grep '^@SQ' | grep -oP 'SN:\\K[^\\s]+'
        )
        if [[ \$bam_ref_name != \$ref_name ]]; then
            >&2 echo -n "BAM file appears to have been aligned against a "
            >&2 echo "different reference:"
            >&2 echo "Ref. name in BAM file: '\$bam_ref_name'"
            >&2 echo "Ref. name in '$reference': '\$ref_name'"
            exit 1
        fi
        # no need to align, just rename
        mv input ${sample_id}.bam
        mv index ${sample_id}.bam.bai
    fi
    stats_from_bam -o ${sample_id}.bamstats -s ${sample_id}.bam.summary -t 4 ${sample_id}.bam
    """
}

process coverageCalc {
    label "wfmpx"
    cpus 1
    memory '2GB'
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
    label "medaka"
    cpus 1
    memory '8GB'
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), val(basecall_model)
        path reference
        path genbank
    output:
        tuple val(sample_id), val(type), path("${sample_id}.annotate.filtered.vcf")
    script:
    """
    medaka consensus ${sample_id}.bam ${sample_id}.hdf --model ${basecall_model}:consensus
    medaka variant --gvcf ${reference} ${sample_id}.hdf ${sample_id}.vcf --verbose
    medaka tools annotate --debug --pad 25 ${sample_id}.vcf ${reference} ${sample_id}.bam ${sample_id}.annotate.vcf
    bcftools filter -e "ALT='.'" ${sample_id}.annotate.vcf | bcftools filter -o ${sample_id}.annotate.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    # vcf-annotator ${sample_id}.annotate.filtered.vcf ${genbank} > ${sample_id}.vcf-annotator.vcf
    """

}

process makeConsensus {
    label "wfmpx"
    cpus 1
    memory '500MB'
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
    cpus 4
    memory '28GB'
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path reference
    output:
        tuple val(sample_id), val(type), path("flye/assembly.fasta"), path("${sample_id}_restricted.fastq"), path(reference)
    script:
    """
    samtools bam2fq ${sample_id}.bam > ${sample_id}_restricted.fastq
    flye --nano-raw ${sample_id}_restricted.fastq  -g 197k -t 4 --meta -o flye
    """
}

process noAssembly{
    label "wfmpx"
    cpus 1
    memory '10MB'
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path reference
    output:
        tuple val(sample_id), val(type), path("medaka/consensus.fasta"), path("${sample_id}_assembly_mapped.bam"), path("${sample_id}_assembly.bed")
    script:
    """
    mkdir medaka
    touch medaka/consensus.fasta
    touch ${sample_id}_assembly_mapped.bam
    touch ${sample_id}_assembly.bed
    """
}

process medaka_polish {
    label "medaka"
    cpus 2
    memory '16GB'
    input:
        tuple val(sample_id),
        val(type),
        path("flye/assembly.fasta"),
        path("${sample_id}_restricted.fastq"),
        path(reference),
        val(basecall_model)
    output:
        tuple val(sample_id), val(type), path("medaka/consensus.fasta"), path("${sample_id}_assembly_mapped.bam")
    script:
    """
    medaka_consensus -i ${sample_id}_restricted.fastq -t 1 -d flye/assembly.fasta \
        -m ${basecall_model}:consensus

    minimap2 -ax map-ont ${reference} medaka/consensus.fasta -t 1 \\
        --cap-kalloc 100m --cap-sw-mem 50m \\
    | samtools sort -o ${sample_id}_assembly_mapped.bam
    """
}

process bamtobed {
    label "wfmpx"
    cpus 1
    memory '500MB'
    input:
        tuple val(sample_id), val(type), path("medaka/consensus.fasta"), path("${sample_id}_assembly_mapped.bam")
    output:
        tuple val(sample_id), val(type), path("medaka/consensus.fasta"), path("${sample_id}_assembly_mapped.bam"), path("${sample_id}_assembly.bed")
    script:
    """
    bedtools bamtobed -i ${sample_id}_assembly_mapped.bam > ${sample_id}_assembly.bed

    """
}

process getVersions {
    label "wfmpx"
    cpus 1
    memory '500MB'
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    flye --version | sed 's/^/flye,/' >> versions.txt
    bcftools --version | head -1 | sed 's/ /,/' >> versions.txt
    """
}

process getVersions_medaka {
    label "medaka"
    cpus 1
    memory '500MB'
    input:
        path "versions.txt"
    output:
        path "versions.txt"
    script:
    """
    medaka --version | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "wfmpx"
    cpus 1
    memory '500MB'
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
    cpus 1
    memory '1GB'
    input:
        path "per-read-stats.tsv.gz"
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
    [ -e ${assembly_bed} ] || echo "" > ${assembly_bed}

    workflow-glue report $report_name \
        --versions versions \
        per-read-stats.tsv.gz \
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
    cpus 1
    memory '500MB'
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
        // get basecall models: we use `params.override_basecaller_cfg` if present;
        // otherwise use `meta.basecall_models[0]` (there should only be one value in
        // the list because we're running ingress with `allow_multiple_basecall_models:
        // false`; note that `[0]` on an empty list returns `null`)
        basecall_models = reads
        | map { meta, reads, index, stats ->
            String basecall_model = \
                params.override_basecaller_cfg ?: meta.basecall_models[0]
            if (!basecall_model) {
                error "Found no basecall model information in the input data for " + \
                    "sample '$meta.alias'. Please provide it with the " + \
                    "`--override_basecaller_cfg` parameter."
            }
            [meta.alias, basecall_model]
        }

        alignment = alignReads(
            reads.map { meta, reads, index, stats -> [meta, reads, index] },
            reference,
        )
        coverage = coverageCalc(alignment.alignment, reference)
        variants = medakaVariants(
            alignment.alignment.join(basecall_models),
            reference,
            genbank
        )
        draft = makeConsensus(variants, reference, coverage)
        software_versions = getVersions()|getVersions_medaka
        workflow_params = getParams()

        if (params.assembly) {
            assembly = flyeAssembly(alignment.alignment, reference)
            | join(basecall_models)
            | medaka_polish
            | bamtobed
        } else {
            assembly = noAssembly(alignment.alignment, reference)
        }

        per_read_stats = reads.map { meta, reads, index, stats ->
            // the glob `"*read*stats.tsv.gz"` matches bamstats and fastcat per-read stats
            [meta, file(stats.resolve("*read*stats.tsv.gz"))]
        }

        report = makeReport(
            per_read_stats.map { meta, stats -> stats },
            software_versions.collect(),
            workflow_params,
            variants.map{it[2]}.collect(),
            coverage,
            reference,
            assembly
        )
    emit:
        per_read_stats
        results = report.concat(
            alignment.alignment.map{it[2]}.collect(),
            alignment.alignment.map{it[3]}.collect(),
            variants.map{it[2]}.collect(),
            variants.map{it[3]}.collect(),
            draft.map{it[2]}.collect(),
            coverage,
            assembly.map{it[2]}.collect(),
            software_versions,
            workflow_params
        )
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    // warn the user if overriding the basecall models found in the inputs
    if (params.override_basecaller_cfg) {
        log.warn \
            "Overriding basecall model with '${params.override_basecaller_cfg}'."
    }

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

    Map ingress_args = [
        "sample": params.sample,
        "sample_sheet": null,
        "stats": true,
        "per_read_stats": true,
        "allow_multiple_basecall_models": false,
    ]

    if (params.fastq) {
        samples = fastq_ingress(ingress_args + [
            "input":params.fastq,
        ])
        | map { meta, reads, stats -> [meta, reads, OPTIONAL_FILE, stats] }
    } else {
        samples = xam_ingress(ingress_args + [
            "input":params.bam,
            "keep_unaligned":true,
        ])
    }

    // make sure we only got one sample
    samples
    | toList
    | subscribe {
        if (it.size() > 1) {
            error "Found more than one barcode / sample; this is currently not supported."
        }
    }

    pipeline(samples, params._reference, params._genbank)

    Channel.empty()
    | mix(
        pipeline.out.per_read_stats
        | map { meta, stats -> [stats, "${meta.alias}.per-read-stats.tsv.gz"] },
        pipeline.out.results
        | map { it -> [it, null] },
    )
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
