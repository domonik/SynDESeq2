

include: "rules/goEnrich.smk"
include: "rules/plots.smk"


rule all:
    input:
        deseqR = expand(rules.GOEnrichment.output, zip, condition=config["conditions"], baseline=config["baselines"]),
        semsim = expand(rules.ClusteredEnrichmentPlot.output, zip, condition=config["conditions"], baseline=config["baselines"], updown=["up", "down"]),
        volcano = expand(rules.volcanoPlot.output, zip, condition=config["conditions"], baseline=config["baselines"]),
        enrich = expand(rules.EnrichmentPlot.output, zip, condition=config["conditions"], baseline=config["baselines"], updown=["up", "down"]),

