from typing import Dict

import numpy as np
import pandas as pd
from plotly import graph_objs as go, express as px

try:
    from example_layout import LAYOUT
except ModuleNotFoundError:
    LAYOUT = go.Layout(
        template="plotly_white",
        font=dict(
            color="black"
        ),
        legend=dict(
            font=dict(
                size=10
            )
        ),
        xaxis=dict(ticklen=0,linecolor="black"),
        yaxis=dict(ticklen=0,linecolor="black"),
        margin=dict(
            b=60,
            t=60,
        )
    )


include: "deseq.smk"
include: "enrich.smk"
import os





rule volcanoPlot:
    input:
        deseq_result = rules.extractDESeqResult.output.result_table,
        uniprot_table = rules.downloadOrganismGOTerms.output.go_terms
    output:
        html=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Volcano/VolcanoPlot_c{condition}_vs_b{baseline}.html"
        ),
        svg=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Volcano/VolcanoPlot_c{condition}_vs_b{baseline}.svg"
        ),
        json=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Volcano/VolcanoPlot_c{condition}_vs_b{baseline}.json"
        ),
    conda: "../envs/DEplots.yml"
    script:
        "../PyScripts/volcanoPlot.py"



def enrichment_plot_from_cp_table(df, mode="scatter"):
    if len(df) == 0:
        return empty_figure()

    def df_div(l):
        return int(l[0]) / int(l[1])
    df["GeneRatio"] = df["GeneRatio"].str.split("/").apply(df_div)
    if mode == "scatter":
        fig = px.scatter(
            df,
            x="GeneRatio",
            y="Description",
            symbol="ONTOLOGY",
            color="p.adjust",
            template="plotly_white",

        )
        fig.update_traces(marker=dict(size=15))

    elif mode == "bar":
        fig = px.bar(
            df,
            x="GeneRatio",
            y="Description",
            color="p.adjust",
            template="plotly_white",
        )
    else:
        raise ValueError(f"mode: {mode} is not valid")
    fig.update_layout(LAYOUT)
    fig.update_layout(
        coloraxis_colorbar=dict(
            yanchor="top",
            y=0.7,
            len=0.7,
            x=1,
            ticks="outside"
        ),
        legend=dict(x=1),
        yaxis=dict(tickmode="linear", type="category", dtick=1)
    )
    return fig


rule EnrichmentPlot:
    input:
        file = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        html = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg"),
        json = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.json")
    run:
        df = pd.read_csv(input.file, sep="\t")
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)
        fig.write_json(output.json)


def reduce_cluster(df: pd.DataFrame, method: str = "strlen", ascending: bool = True):
    non_clustered = df[df.Cluster == -1]
    df = df[~(df.Cluster == -1)]
    if method == "strlen":
        df["sort_column"] = df["Description"].str.len()
        ascending = ascending
    elif method in df.columns:
        df["sort_column"] = df[method]
        ascending = ascending
    else:
        raise NotImplementedError(f"Method '{method}' not implemented")
    df = df.sort_values("sort_column", ascending=ascending)
    df = df.drop_duplicates(["ONTOLOGY", "Cluster"], keep='first')
    df = pd.concat((df, non_clustered), axis=0)
    df = df.drop("sort_column", axis=1)
    return df


rule ClusteredEnrichmentPlot:
    input:
        file = rules.ClusterSemSim.output.table,
    output:
        html=os.path.join(config[
            "RUN_DIR"],"PipelineData/Plots/Enrichment/ClusteredGOEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Enrichment/ClusteredGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg"),
        json=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Enrichment/ClusteredGOEnrichment_{updown}_c{condition}_vs_b{baseline}.json")
    run:
        df = pd.read_csv(input.file, sep="\t")
        df = reduce_cluster(df, method=config["sort_method"], ascending=config["ascending"])
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)
        fig.write_json(output.json)

rule KEGGEnrichmentPlot:
    input:
        file = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        html = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.svg"),
        json = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.json")
    run:
        df = pd.read_csv(input.file, sep="\t")
        df["ONTOLOGY"] = "KEGG"
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)
        fig.write_json(output.json)


def empty_figure():
    fig = go.Figure()
    fig.add_annotation(
        xref="paper",
        yref="paper",
        xanchor="center",
        yanchor="middle",
        x=0.5,
        y=0.5,
        text="No Data to display",
        showarrow=False,
        font=(dict(size=16))
    )
    return fig


def add_boxes(fig: go.Figure, config):
    fig.add_shape(
        type="rect",
        x0=config["log2FCCutOff"], y0=-1 * np.log10(config["pAdjCutOff"]), x1=100, y1=200,
        line=dict(
            color="grey",
            width=2,
        ),
        layer="below",
        fillcolor=config["upColor"],
        opacity=0.05
    )
    fig.add_shape(
        type="rect",
        x0=-config["log2FCCutOff"], y0=-1 * np.log10(config["pAdjCutOff"]), x1=-100, y1=200,
        line=dict(
            color="grey",
            width=2,
        ),
        layer="below",
        fillcolor=config["downColor"],
        opacity=0.05
    )
    fig.add_hline(y=-1 * np.log10(config["pAdjCutOff"]), line_color="grey", layer="below", line_dash="dash")
    fig.add_vline(x=-config["log2FCCutOff"], line_color="grey", layer="below", line_dash="dash")
    fig.add_vline(x=config["log2FCCutOff"], line_color="grey", layer="below", line_dash="dash")
    return fig


def compare_results(files, sep, condition, base, padj_cutoff, lfc_cutoff, colors = ("rgb(0, 142, 151) ", "rgb(252, 76, 2)", "#2CA02C"), y: str = "Name"):
    fig = go.Figure()
    df = {"Enriched in": [], "Percentage": [], y: []}
    for name, file in files:
        df1 = pd.read_csv(file, sep=sep)
        df1 = df1.dropna()
        df1_up = df1[(df1["padj"] <= padj_cutoff) & (df1["log2FoldChange"] >= lfc_cutoff)]
        df1_down = df1[(df1["padj"] <= padj_cutoff) & (df1["log2FoldChange"] <= -lfc_cutoff)]
        lup = len(df1_up)
        ldown = len(df1_down)
        lges = len(df1)

        df["Percentage"] += [(ldown/lges) * 100, ((lges - lup -ldown) / lges) * 100, (lup/lges) * 100]
        df[y] += [name, name, name]
        df["Enriched in"] += [base, "not Enriched", condition]
    df = pd.DataFrame(df)
    fig = px.bar(df, x="Percentage", y=y, color="Enriched in", color_discrete_sequence=colors)
    fig.update_traces(marker_line_color='black',
                      marker_line_width=1)
    fig.update_xaxes(range=[0, 100])

    return fig
