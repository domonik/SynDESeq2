import os

include: "deseq.smk"


ORGDB = os.path.join(config["GOTermDir"], "AnnotationDBs", ORGANISMID)

rule clusterProfilerInstallFromGitHub:
    # This is necessary since they don´t update their conda package
    output:
        lib = directory(os.path.join(config["Rlib"], "clusterProfiler"))
    conda:
        "../envs/REnvironment.yml"
    script:
        "../Rscripts/installClusterProfiler.R"


rule buildLocalKEGGdb:
    input: rules.clusterProfilerInstallFromGitHub.output
    params:
        localKEGGdb = os.path.join(config["Rlib"], config["keggOrgID"])
    output:
        lib = os.path.join( config["Rlib"] , config["keggOrgID"], "KEGG.db_1.0.tar.gz")
    conda:
        "../envs/REnvironment.yml"
    script:
        "../Rscripts/downloadKEGGdb.R"

rule GOEnrichment:
    input:
        #cp = rules.clusterProfilerInstallFromGitHub.output.lib,
        annotation_db = ORGDB,
        deseq_results = rules.extractDESeqResult.output.result_table
    conda:
        "../envs/REnvironment.yml"
    output:
        up = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_up_c{condition}_vs_b{baseline}.tsv"),
        down = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_down_c{condition}_vs_b{baseline}.tsv")
    script:
        "../Rscripts/goEnrichment.R"


rule GSEAGO:
    input:
        #cp = rules.clusterProfilerInstallFromGitHub.output.lib,
        annotation_db = ORGDB,
        deseq_results = rules.extractDESeqResult.output.result_table,
    conda:
        "../envs/REnvironment.yml"
    output:
        enriched = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GSEAGO_c{condition}_vs_b{baseline}.tsv"),
        gsdata = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GSEAGO_plot_data_c{condition}_vs_b{baseline}.tsv"),
    script:
        "../Rscripts/GSEA.R"


rule SemanticSimilarity:
    input:
        annotation_db = ORGDB,
        enrichData = os.path.join(
            config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    output:
        table = os.path.join(
            config["RUN_DIR"],"PipelineData/Enrichment/SemanticSimilarity_ont{subcat}_{updown}_c{condition}_vs_b{baseline}.tsv"
        )
    conda:
        "../envs/REnvironment.yml"
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
                "RUN_DIR"],"PipelineData/Enrichment/ClusteredGOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv"
            )
    conda: "../envs/sklearn.yml"
    script:  "../PyScripts/clusterGOEnrich.py"




rule enrichKEGG:
    input:
        #cp = rules.clusterProfilerInstallFromGitHub.output.lib,
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
    conda:
        "../envs/REnvironment.yml"
    script:
        "../Rscripts/keggEnrichment.R"
