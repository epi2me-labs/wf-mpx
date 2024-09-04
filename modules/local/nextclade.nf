process nextclade {
    label "nextclade"
    cpus 1
    input:
        path consensus
        val nextclade_type
    output:
        path("nextclade_${nextclade_type}.json"), optional: true
    script:
    """
    nextclade run --dataset-name MPXV --output-all=nextclade_output/ ${consensus}
    mv nextclade_output/nextclade.json nextclade_${nextclade_type}.json
    """
}
