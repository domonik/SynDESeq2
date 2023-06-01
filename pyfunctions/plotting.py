import pandas as pd
import plotly.graph_objs as go
import numpy as np
from typing import Dict
from example_layout import LAYOUT
import plotly.express as px


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
            text=df["plot_name"],
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


def enrichment_plot_from_cp_table(df):
    if len(df) == 0:
        return empty_figure()

    def df_div(l):
        return int(l[0]) / int(l[1])
    df["GeneRatio"] = df["GeneRatio"].str.split("/").apply(df_div)

    fig = px.scatter(
        df,
        x="GeneRatio",
        y="Description",
        symbol="ONTOLOGY",
        color="p.adjust",
        template="plotly_white",

    )
    fig.update_traces(marker=dict(size=15))
    fig.update_layout(
        coloraxis_colorbar=dict(
            yanchor="top",
            y=0.7,
            len=0.7,
            x=1,
            ticks="outside"
        ),
        legend=dict(x=1),
        yaxis=dict(tickmode="linear")
    )
    fig.update_layout(LAYOUT)
    return fig

