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




