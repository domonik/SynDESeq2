import os
from pyfunctions.helpers import cluster_go_enrich
include: "setup.smk"
include: "deseq.smk"




rule clusterProfilerInstallFromGitHub:
    # This is necessary since they don´t update their conda package
    output:
        lib = directory(os.path.join(config["Rlib"], "clusterProfiler"))
    script:
        "../Rscripts/installClusterProfiler.R"


rule buildLocalKEGGdb:
    input: rules.clusterProfilerInstallFromGitHub.output
    params:
        localKEGGdb = os.path.join(config["Rlib"], config["keggOrgID"])
    output:
        lib = os.path.join( config["Rlib"] , config["keggOrgID"], "KEGG.db_1.0.tar.gz")
    script:
        "../Rscripts/downloadKEGGdb.R"

rule GOEnrichment:
    input:
        cp = rules.clusterProfilerInstallFromGitHub.output.lib,
        annotation_db = rules.generateOrgDB.output.annotation_db,
        deseq_results = rules.extractDESeqResult.output.result_table
    output:
        up = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_up_c{condition}_vs_b{baseline}.tsv"),
        down = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_down_c{condition}_vs_b{baseline}.tsv")
    script:
        "../Rscripts/goEnrichment.R"


rule SemanticSimilarity:
    input:
        annotation_db = rules.generateOrgDB.output.annotation_db,
        enrichData = os.path.join(
            config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    output:
        table = os.path.join(
            config["RUN_DIR"],"PipelineData/Enrichment/SemanticSimilarity_ont{subcat}_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    script:
        "../Rscripts/calcSemSim.R"


rule ClusterSemSim:
    input:
        semsim = expand(
            rules.SemanticSimilarity.output.table, subcat=["MF", "BP", "CC"], allow_missing=True
        ),
        enrichment = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        table = os.path.join(
            config[
                "RUN_DIR"],"PipelineData/Enrichment/ClusteredEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
            )
    run:
        df = cluster_go_enrich(input.semsim, input.enrichment, config=config)
        df.to_csv(output.table, sep="\t", index=False)



rule enrichKEGG:
    input:
        cp = rules.clusterProfilerInstallFromGitHub.output.lib,
        defile = rules.extractDESeqResult.output.result_table,
    output:
        up = os.path.join(
            config[
                "RUN_DIR"],"PipelineData/Enrichment/KEGGEnrichment_up_c{condition}_vs_b{baseline}.tsv"
            ),
        down = os.path.join(
            config[
                "RUN_DIR"],"PipelineData/Enrichment/KEGGEnrichment_down_c{condition}_vs_b{baseline}.tsv"
            )
    script:
        "../Rscripts/keggEnrichment.R"
