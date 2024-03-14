from typing import Dict

import numpy as np
import pandas as pd
from plotly import graph_objs as go, express as px

from example_layout import LAYOUT


include: "deseq.smk"
include: "enrich.smk"
import os


def volcano_from_deseq_result(
        deseq_result,
        config: Dict,
        initial_sep: str = ",",
        tag2name: str = None,
        highlight: Dict = None
):
    hovertemplate = '<i>Y</i>: %{y:.2f}' + \
                    '<br><b>X</b>: %{x}<br>' + \
                    '<b>%{text}</b>'
    df = pd.read_csv(deseq_result, sep="\t", index_col=0)
    if tag2name:
        tag2name = pd.read_csv(tag2name, sep=initial_sep, index_col=0)
        df = pd.concat((df, tag2name), axis=1)
        df["plot_name"] = df["gene_name"].str.cat(df.index, sep="-")
        df = df.dropna()
    else:
        df["plot_name"] = df.index
        df["gene_name"] = df.index
    df["-log10padj"] = -1 * np.log10(df["padj"])
    max_log10padj = np.ceil(df["-log10padj"].max())
    min_fc = np.floor(df["log2FoldChange"].min())
    max_fc = np.ceil(df["log2FoldChange"].max())

    fig = go.Figure(layout=LAYOUT)

    if highlight is not None:
        for key, value in highlight.items():
            color, names = value
            if isinstance(names, str):
                mask = df.gene_name.str.contains(names)
                to_highlight = df[mask]
                n = ~mask
            else:
                to_highlight = df[df.index.isin(names)]
                if len(to_highlight) == 0:
                    to_highlight = df[df.gene_name.isin(names)]
                    n = ~df["gene_name"].isin(names)
                else:
                    n = ~df.index.isin(names)
            df = df[n]
            fig.add_trace(go.Scatter(
                x=to_highlight["log2FoldChange"],
                y=to_highlight["-log10padj"],
                mode="markers",
                marker=dict(color=color),
                hovertemplate=hovertemplate,
                text=to_highlight["plot_name"],
                name=key
            ))
    df["significant"] = "not"
    df.loc[(df["log2FoldChange"] >= config["log2FCCutOff"]) & (df["padj"] < config["pAdjCutOff"]), "significant"] = "up"
    df.loc[(df["log2FoldChange"] <= -config["log2FCCutOff"]) & (df["padj"] < config["pAdjCutOff"]), "significant"] = "down"

    for sig, color in {"up": config["upColor"], "down": config["downColor"], "not": config["normalColor"]}.items():

        fig.add_trace(go.Scatter(
            x=df[df["significant"] == sig]["log2FoldChange"],
            y=df[df["significant"] == sig]["-log10padj"],
            mode="markers",
            marker=dict(color=color),
            hovertemplate=hovertemplate,
            text=df[df["significant"] == sig]["plot_name"],
            name="Normal Genes",
            showlegend=False
        ))
    fig.data = fig.data[::-1]
    fig = add_boxes(fig, config)
    fig.update_xaxes(
        range=[min_fc, max_fc],
        dtick=1
    )
    fig.update_yaxes(
        range=[-0.5, max_log10padj],
    )
    fig.update_layout(
        xaxis_title="Log2FoldChange",
        yaxis_title="-Log10(pval)"
    )
    return fig


rule volcanoPlot:
    input:
        deseq_result = rules.extractDESeqResult.output.result_table
    output:
        html=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Volcano/VolcanoPlot_c{condition}_vs_b{baseline}.html"
        ),
        svg=os.path.join(
            config["RUN_DIR"],"PipelineData/Plots/Volcano/VolcanoPlot_c{condition}_vs_b{baseline}.svg"
        ),
    run:
        fig = volcano_from_deseq_result(
            deseq_result=input.deseq_result,
            initial_sep=config["initial_sep"],
            config=config,
            tag2name=config["tag2Name"],
            highlight=config["highlight"]
        )
        title = config["design"] + f" - {wildcards.baseline} vs {wildcards.condition}"
        fig.update_layout(title=dict(text=title))
        fig.write_html(output.html)
        fig.write_image(output.svg)


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
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg")
    run:
        df = pd.read_csv(input.file, sep="\t")
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)


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
            config["RUN_DIR"],"PipelineData/Plots/Enrichment/ClusteredGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg")
    run:
        df = pd.read_csv(input.file, sep="\t")
        df = reduce_cluster(df, method=config["sort_method"], ascending=config["ascending"])
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)

rule KEGGEnrichmentPlot:
    input:
        file = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        html = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/KEGGEnrichment_{updown}_c{condition}_vs_b{baseline}.svg")
    run:
        df = pd.read_csv(input.file, sep="\t")
        df["ONTOLOGY"] = "KEGG"
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)


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
