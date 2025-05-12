import os
import yaml
from snakemake.utils import min_version
min_version("6.0")

configfile: "multilevelConfig.yml"

with open(config["primary-file"], "r") as handle:
    config["primary"] = yaml.safe_load(handle)

with open(config["secondary-file"], "r") as handle:
    config["secondary"] = yaml.safe_load(handle)


module other_workflow:
    snakefile: "snakefile.smk"
    config: config["primary"]

use rule * from other_workflow as first_*


rule runOnLeastDiffRegGenes:
    input:
        files = expand(rules.first_extractDESeqResult.output.result_table,  zip, condition=config["primary"]["conditions"], baseline=config["primary"]["baselines"])
    output:
        file =   os.path.join(config["secondary"]["RUN_DIR"], "PipelineData/Multilevel/intersectingGenes.tsv")
    run:
        import pandas as pd
        sets = []
        for file in input.files:
            df = pd.read_csv(file, sep="\t")
            df = df[(df["log2FoldChange"] <= config["multilevel"]["lfc_cutoff"]) & (df["log2FoldChange"] >= -config["multilevel"]["lfc_cutoff"]) ]
            df = df[df["padj"] >= config["multilevel"]["padj_cutoff"]]
            sets.append(set(df.index))
        hs_genes = config["multilevel"]["additional_housekeeping_genes"]
        if len(hs_genes) > 0:
            for gene in hs_genes:
                assert gene in df.index, f"{gene} not found in index"
            sets.append(hs_genes)
        if config["multilevel"]["mode"] == "intersect":
            sets = set.intersection(*sets)
        else:
            raise NotImplementedError("other mode not implemented yet")
        df = df[df.index.isin(sets)]
        df.to_csv(output.file, sep="\t")

config["secondary"]["housekeeping"] = str(rules.runOnLeastDiffRegGenes.output.file)
config["secondary"]["use-housekeeping"] = True

module another_workflow:
    snakefile: "snakefile.smk"
    config: config["secondary"]



use rule * from another_workflow as second_*


rule plotOverlap:
    input:
        deseq1 = rules.first_extractDESeqResult.output.result_table,
        deseq2 = rules.second_extractDESeqResult.output.result_table
    output:
        tsv = os.path.join(config["secondary"]["RUN_DIR"], "PipelineData/Multilevel/c{condition}_vs_b{baseline}_overlap.tsv")
    run:
        import pandas as pd
        des1 = pd.read_csv(input.deseq1, sep="\t")
        des2 = pd.read_csv(input.deseq2, sep="\t")
        des1 = des1[(des1["log2FoldChange"] >= config["primary"]["log2FCCutOff"]) & (des1["padj"] >= config["primary"]["pAdjCutOff"])]

        des2 = des2[(des2["log2FoldChange"] >= config["primary"]["log2FCCutOff"]) & (des2["padj"] >= config["primary"]["pAdjCutOff"])]
        print(des1)
        print("_________________")
        print(des2)
        des = des1.join(des2, lsuffix="_first", rsuffix="_second", how="inner")
        print(des)

rule all:
    default_target: True
    input:
        first=rules.first_all.input,
        second=rules.second_all.input,
        overlap = expand(rules.plotOverlap.output, zip, condition=config["primary"]["conditions"], baseline=config["primary"]["baselines"])


