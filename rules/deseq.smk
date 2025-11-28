import os

include: "setup.smk"

rule extractDESeqResult:
    input:
        result = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DESeqResult.RData"),
    output:
        result_table = os.path.join(config["RUN_DIR"], "PipelineData/DESeqResults/DESeqResult_c{condition}_vs_b{baseline}.tsv")
    params:
        factor = config["design"].split("+")[-1].split(" ")[-1]
    conda:
        "../envs/DESeq2.yml"
    script:
        "../Rscripts/extractDESeqResult.R"


rule plotNormFactors:
    input:
        file = rules.dropUnusedSamples.output.counts
    output:
        file = os.path.join(config["RUN_DIR"], "PipelineData/Plots/DESeq/NormFactorsHisto.png")
    run:
        import pandas as pd
        import numpy as np
        import plotly.graph_objs as go
        from plotly.colors import DEFAULT_PLOTLY_COLORS
        from plotly.subplots import make_subplots
        df = pd.read_csv(input.file, sep="\t", index_col=0)
        n = df.shape[1]
        df = df.astype(float)
        x = df.product(axis=1)
        df["pseudoref"] = np.power(df.product(axis=1), (1./n))
        for col in df.columns:
            if col != 'pseudoref':  # Skip the fixed column itself
                new_col_name = col + '_ratio'
                df[new_col_name] = df[col] / df["pseudoref"]
        df = df.replace([np.inf, -np.inf], np.nan)
        cols = [col for col in df.columns if "_ratio" in col]
        fig = go.Figure()
        for idx, col in enumerate(cols):
            array = df[col]
            fig.add_trace(
                go.Histogram(
                    x=array,
                    name=col,
                )
            )
            fig.add_vline(x=np.nanmedian(array), line=dict(color=DEFAULT_PLOTLY_COLORS[idx]))
        fig.update_layout(barmode='group',bargap=0,bargroupgap=0.0)
        fig.show()
        raise

rule runDESeq:
    input:
        counts = rules.dropUnusedSamples.output.counts,
        annotation = rules.dropUnusedSamples.output.annotation,
        housekeeping = branch(config["use-housekeeping"], then=config["housekeeping"], otherwise=None)
    params:
        use_housekeeping =  config["use-housekeeping"],
        housekeeping= config["housekeeping"],
        use_spike_ins = config["use-spike-ins"],
        spike_in_pattern = config["spike-in-pattern"],
        design = config["design"]
    conda:
        "../envs/DESeq2.yml"
    output:
        heatmap = os.path.join(config["RUN_DIR"], "PipelineData/Plots/DESeq/SampleHeatmap.png"),
        pca_data = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/PCAData.tsv"),
        corr_data = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/CorrData.tsv"),
        dispersion_estimates = os.path.join(config["RUN_DIR"], "PipelineData/Plots/DESeq/DispersionEstimates.png"),
        deseq_result = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DESeqResult.RData"),
        normalized_counts = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/normalized_counts.tsv"),
        size_factors = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/size_factors.tsv")
    script:
        "../Rscripts/deseq.R"




