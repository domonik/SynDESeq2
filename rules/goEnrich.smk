import os
from pyfunctions.helpers import cluster_go_enrich
include: "setup.smk"
include: "deseq.smk"


rule GOEnrichment:
    input:
        annotation_db = rules.generateOrgDB.output.annotation_db,
        deseq_results = rules.extractDESeqResult.output.result_table
    output:
        up = os.path.join(config["RUN_DIR"], "PipelineData/GOEnrichment/GOEnrichment_up_c{condition}_vs_b{baseline}.tsv"),
        down = os.path.join(config["RUN_DIR"], "PipelineData/GOEnrichment/GOEnrichment_down_c{condition}_vs_b{baseline}.tsv")
    script:
        "../Rscripts/goEnrichment.R"


rule SemanticSimilarity:
    input:
        annotation_db = rules.generateOrgDB.output.annotation_db,
        enrichData = os.path.join(
            config["RUN_DIR"], "PipelineData/GOEnrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    output:
        table = os.path.join(
            config["RUN_DIR"],"PipelineData/GOEnrichment/SemanticSimilarity_ont{subcat}_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    script:
        "../Rscripts/calcSemSim.R"


rule ClusterSemSim:
    input:
        semsim = expand(
            rules.SemanticSimilarity.output.table, subcat=["MF", "BP", "CC"], allow_missing=True
        ),
        enrichment = os.path.join(config["RUN_DIR"], "PipelineData/GOEnrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        table = os.path.join(
            config[
                "RUN_DIR"],"PipelineData/GOEnrichment/ClusteredEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
            )
    run:
        df = cluster_go_enrich(input.semsim, input.enrichment, config=config)
        df.to_csv(output.table, sep="\t", index=False)

