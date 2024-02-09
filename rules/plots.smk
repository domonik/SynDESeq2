import pandas as pd


include: "deseq.smk"
include: "enrich.smk"
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
        from pyfunctions.plotting import volcano_from_deseq_result

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
        file = os.path.join(config["RUN_DIR"], "PipelineData/Enrichment/GOEnrichment_{updown}_c{condition}_vs_b{baseline}.tsv")
    output:
        html = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "PipelineData/Plots/Enrichment/RawGOEnrichment_{updown}_c{condition}_vs_b{baseline}.svg")
    run:
        from pyfunctions.plotting import enrichment_plot_from_cp_table

        df = pd.read_csv(input.file, sep="\t")
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
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
        from pyfunctions.plotting import enrichment_plot_from_cp_table
        from pyfunctions.helpers import reduce_cluster


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
        from pyfunctions.plotting import enrichment_plot_from_cp_table

        df = pd.read_csv(input.file, sep="\t")
        df["ONTOLOGY"] = "KEGG"
        fig = enrichment_plot_from_cp_table(df, mode=config["enrichPlotType"])
        fig.write_html(output.html)
        fig.write_image(output.svg)