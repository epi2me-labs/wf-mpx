#!/usr/bin/env python
"""Create workflow report."""
import json
import math
import os

from aplanat import bars
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import ont_colors
from bokeh.models import BasicTickFormatter, ColumnDataSource, LabelSet
from bokeh.plotting import figure
import pandas as pd
import pysam
import vcf

from .util import wf_parser  # noqa: ABS101


def load_fasta(reference):
    """Load reference data."""
    fasta = pysam.Fastafile(reference)
    return fasta


def load_vcf(vcf_file):
    """Load VCF data."""
    vcf_reader = vcf.Reader(filename=vcf_file)

    # define our data structure
    data = dict(
        chromosome=list(),
        start=list(),
        end=list(),
        quality=list(),
        reference=list(),
        alternate=list(),
        is_indel=list()
    )

    info_fields = ['DP', 'DPS']

    for info_field in info_fields:
        data[info_field] = list()

    # for records in VCF file add to our structure
    for record in vcf_reader:
        for alt in record.ALT:
            data['chromosome'].append(record.CHROM)
            data['start'].append(record.POS)
            data['end'].append(record.POS+len(alt))
            data['quality'].append(record.QUAL)
            ref = record.REF
            alt = str(alt)
            if len(ref) > 10:
                ref = f"{ref[0:10]}...({len(ref)})"
            data['reference'].append(ref)
            if len(alt) > 10:
                alt = f"{alt[0:10]}...({len(alt)})"
            data['alternate'].append(alt)
            data['is_indel'].append(record.is_indel)

            for info_field in info_fields:
                field_data = record.INFO[info_field]
                if isinstance(field_data, list):
                    field_data = ",".join(str(x) for x in field_data)
                data[info_field].append(field_data)

    df = pd.DataFrame.from_dict(data)
    df['width'] = df['end']-df['start']

    return df


def make_coverage_section(coverage_file, variants_data, report_doc):
    """Make the coverage section."""
    coverage = pd.read_csv(coverage_file, sep="\t", header=0)

    section = report_doc.add_section()
    section._add_item("""<h3>Genome coverage</h3>

    <p>The plot below shows read depth across the whole genome, annotated with variants
    in the reference-based assembly called by Medaka. Dots represent single nucleotide
    polymorphisms, while bars represent insertions and deletions.</p>""")

    max_coverage = max(coverage['depth'])

    p = figure(
            x_axis_label='position',
            y_axis_label='depth',
            width=1000,
            height=300)

    # add a line for the depth of coverage
    p.line(
            coverage['pos'],
            coverage['depth'],
            line_color=ont_colors.BRAND_BLUE)

    # if VCF is non-empty
    if not variants_data.empty:
        # add a point for the snps
        p.circle(
            variants_data[~variants_data.is_indel]['start'],
            y=5,
            size=5,
            color=ont_colors.BRAND_GREY)

        # add a bar for the indels
        p.hbar(
            left=variants_data[variants_data.is_indel]['start'],
            y=6,
            right=variants_data[variants_data.is_indel]['end'],
            height=max_coverage/10,
            color=ont_colors.BRAND_LIGHT_BLUE)

    p.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    section.plot(p)


def get_context(variant, fasta):
    """Get the sequence context around the variant of interest."""
    start = variant[0]
    seq = fasta.fetch(
        reference=fasta.references[0], start=start-1, end=start+1)

    if len(variant[1]) > 1:
        return "NA"
    base = fasta.fetch(
        reference=fasta.references[0], start=start, end=start+1)

    return f"""{seq}>{variant[1]}{base}"""


def make_variants_context(variants_data, reference, report):
    """Make variants context section."""
    section = report.add_section()
    section._add_item("""<h3>Mutational spectrum</h3>

    <p>Evidence of APOBEC3 host enzyme editing has been observed in Monkeypox virus
    (MPXV), as described in <a href="https://tinyurl.com/2fawchvu">this virological.org
    post</a> by Áine O’Toole & Andrew Rambaut. The mutational spectrum of the
    reference-based assembly is shown below.</p>""")

    variants_data['context'] = [
        get_context(x, reference) for x in zip(
            variants_data['start'], variants_data['alternate'])]

    summary = variants_data.context.value_counts().to_frame().set_axis(
        ['Count'], axis=1)

    summary.index.name = 'Variant'
    summary = summary.reset_index()
    summary = summary[summary.Variant != "NA"]

    variants_context_plot = bars.simple_bar(
        summary['Variant'], summary['Count'])
    variants_context_plot.xaxis.major_label_orientation = math.pi/4
    variants_context_plot.xaxis.axis_label = 'Variant'
    variants_context_plot.yaxis.axis_label = 'Count'
    section.plot(variants_context_plot)


def make_variants_table(variants_data, report):
    """Make the variants table."""
    section = report.add_section()
    section._add_item("""<h3>All variants</h3>

    <p>The following table shows variants in the reference-based assembly called
    by Medaka with respect to the chosen reference sequence.</p>""")
    section.table(
            variants_data,
            index=False,
            th_color=ont_colors.BRAND_BLUE,
            paging=True,
            searchable=False)
    pass


def make_assembly_summary(bed, report):
    """Make a plot for the assembly."""
    section = report.add_section()
    section._add_item("""<h3>De novo assembly</h3>

    <p>The plot below shows an alignment of <em>de novo</em> assembly contigs to the
    chosen reference sequence.</p>""")

    contigs = pd.read_csv(bed, sep='\t', header=0).set_axis(
            ['reference', 'start', 'end', 'name', 'qual', 'strand'],
            axis=1)

    p = figure(
            x_axis_label='position',
            y_axis_label='contig',
            x_range=(-1000, 200000),
            width=1000,
            height=300)

    p.hbar(
        left=contigs.start,
        right=contigs.end,
        y=contigs.index,
        height=1,
        color=ont_colors.BRAND_BLUE)
    contigs['text_position'] = (contigs.start-1000)
    labels = LabelSet(
        x='text_position',
        y='index',
        text='name',
        text_align='right',
        text_baseline='middle',
        source=ColumnDataSource(contigs),
        text_color=ont_colors.BRAND_BLUE)
    p.add_layout(labels)

    p.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    section.plot(p)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--variants", required=True,
        help="medaka variants VCF file")
    parser.add_argument(
        "--reference", required=True,
        help="reference fasta file")
    parser.add_argument(
        "--coverage", required=True,
        help="depth of coverage file")
    parser.add_argument(
        "--assembly_bed", required=True,
        help="bed file of assembly mapped to reference")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--nextclade_ref", default=None, required=True,
        help="A JSON file containing the nextclade result")
    parser.add_argument(
        "--nextclade_denovo", default=None,
        help="A JSON file containing the nextclade result")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    return parser


def main(args):
    """Run the entry point."""
    report = WFReport(
        "wf-mpx Sequencing Report", "wf-mpx",
        revision=args.revision, commit=args.commit)

    section = report.add_section()
    section.markdown(""" ### Preamble

    This workflow is intended to assemble a __draft__ consensus Monkeypox virus (MPXV)
    genome from Oxford Nanopore Technologies Sequencing data.

    The rough outline of this workflow:

    * Map reads to a reference genome
    * Call variants with respect to this reference using Medaka
    * Create a consensus sequence by applying these variants to the reference
    * Attempt *de novo* assembly with Flye, polished with Medaka (optional)

    Assembled consensus sequences are masked ('N') in regions of <20x coverage,
    deletions are represented as "-", and insertions are in lower case. Resulting
    assemblies __will require manual review__.""")

    if os.path.exists(args.summaries[0]):
        report.add_section(
            section=fastcat.full_report(args.summaries))

    with open(args.nextclade_ref, 'r') as file:
        data = json.load(file)

    section = report.add_section()
    section.markdown("""
### Nextclade results
#### Reference-based assembly
""")
    if len(data['results']) > 0:
        section.table(pd.DataFrame(
            {
                'clade': [data['results'][0]['clade']],
                'coverage': [data['results'][0]['coverage']],
                'overall_qc': [data['results'][0]['qc']['overallStatus']]}))
    else:
        section.markdown("""
        Clade typing was unsuccessful.
        """)

    if len(data['errors']) > 0:
        section.table(pd.DataFrame(data['errors']))

    if args.nextclade_denovo:
        with open(args.nextclade_denovo, 'r') as file:
            data = json.load(file)

        section.markdown("""
        #### *De novo* assembly
        """)
        if len(data['results']) > 0:
            section.table(pd.DataFrame(
                {
                    'clade': [data['results'][0]['clade']],
                    'coverage': [data['results'][0]['coverage']],
                    'overall_qc': [data['results'][0]['qc']['overallStatus']]}))
        else:
            section.markdown("""
            Clade typing was unsuccessful.
            """)

        if len(data['errors']) > 0:
            section.table(pd.DataFrame(data['errors']))

    variants = load_vcf(args.variants)

    make_coverage_section(args.coverage, variants, report)

    reference = load_fasta(args.reference)

    make_variants_context(variants, reference, report)

    make_variants_table(variants, report)

    if not open(args.assembly_bed).readline() == '':
        make_assembly_summary(args.assembly_bed, report)

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)
