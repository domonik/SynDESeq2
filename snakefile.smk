

include: "rules/enrich.smk"
include: "rules/plots.smk"


rule all:
    input:
        deseqR = expand(rules.GOEnrichment.output, zip, condition=config["conditions"], baseline=config["baselines"]),
        volcano = expand(rules.volcanoPlot.output, zip, condition=config["conditions"], baseline=config["baselines"]),
        enrich = expand(expand(rules.EnrichmentPlot.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True), updown=["up", "down"]),
        enrich2 = expand(expand(rules.ClusteredEnrichmentPlot.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True), updown=["up", "down"]),
        enrichkegg = expand(expand(rules.KEGGEnrichmentPlot.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True), updown=["up", "down"]),

