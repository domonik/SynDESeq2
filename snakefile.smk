
#configfile:  "config.yml"

include: "rules/enrich.smk"



rule all:
    input:
        deseqR = expand(rules.GOEnrichment.output, zip, condition=config["conditions"], baseline=config["baselines"]),
        enrich2 = expand(expand(rules.ClusterSemSim.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True), updown=["up", "down"]),
        enrichkegg = expand(expand(rules.enrichKEGG.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True), updown=["up", "down"]),
        gsea_go = expand(expand(rules.GSEAGO.output, zip, condition=config["conditions"], baseline=config["baselines"], allow_missing=True)),

