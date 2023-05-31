import os

include: "setup.smk"




rule runDESeq2:
    input:
        counts = rules.dropUnusedSamples.output.counts,
        annotation = rules.dropUnusedSamples.output.annotation,
    output:
        heatmap = os.path.join(config["RUN_DIR"], "PipelineData/Plots/DESeq/SampleHeatmap.png"),
        pca_data = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/PCAData.tsv"),
        dispersion_estimates = os.path.join(config["RUN_DIR"], "PipelineData/Plots/DESeq/DispersionEstimates.png"),
        deseq_result = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DESeqResult.RData")
    script:
        "../Rscripts/deseq.R"


rule extractDESeqResult:
    input:
        result = rules.runDESeq2.output.deseq_result
    output:
        result_table = os.path.join(config["RUN_DIR"], "PipelineData/DESeqResults/DESeqResult_c{condition}_vs_b{baseline}.tsv")
    params:
        factor = config["design"].split("+")[-1].split(" ")[-1]
    script:
        "../Rscripts/extractDESeqResult.R"

