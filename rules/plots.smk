import pandas as pd

from pyfunctions.plotting import volcano_from_deseq_result, enrichment_plot_from_cp_table
from pyfunctions.helpers import reduce_cluster

include: "deseq.smk"
include: "goEnrich.smk"
import os



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


rule EnrichmentPlot:
    input:
        file = os.path.join(config["RUN_DIR"], "PipelineData/GOEnrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        html = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg")
    run:
        df = pd.read_csv(input.file, sep="\t")
        fig = enrichment_plot_from_cp_table(df)
        fig.write_html(output.html)
        fig.write_image(output.svg)


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
        fig = enrichment_plot_from_cp_table(df)
        fig.write_html(output.html)
        fig.write_image(output.svg)
