import os

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

include: "setup.smk"
include: "deseq.smk"




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
        cp = rules.clusterProfilerInstallFromGitHub.output.lib,
        annotation_db = rules.generateOrgDB.output.annotation_db,
        deseq_results = rules.extractDESeqResult.output.result_table
    conda:
        "../envs/REnvironment.yml"
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
    conda:
        "../envs/REnvironment.yml"
    script:
        "../Rscripts/calcSemSim.R"


def cluster_go_enrich(files, enrich_file, config):
    enrich_df = pd.read_csv(enrich_file, sep="\t")
    if len(enrich_df) == 0:
        enrich_df["Cluster"] = None
        return enrich_df

    rows = []
    clusters = []
    for file in files:
        df = pd.read_csv(file, sep="\t", index_col=0)
        if np.isnan(df.iloc[0, 0]):
            continue
        distance_matrix = 1 - np.asarray(df.iloc[:, :])
        clustering = DBSCAN(metric="precomputed", eps=config["cluster-eps"], min_samples=2)
        fitted = clustering.fit_predict(distance_matrix)
        rows += list(df.columns)
        clusters += list(fitted)

    clusters = pd.DataFrame(
        {
            "ID": rows,
            "Cluster": clusters
        }

    )
    df = pd.merge(enrich_df, clusters, on="ID")
    return df


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
    conda:
        "../envs/REnvironment.yml"
    script:
        "../Rscripts/keggEnrichment.R"
