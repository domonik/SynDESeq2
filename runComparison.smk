import os
import yaml
from snakemake.utils import min_version
min_version("6.0")

configfile: "runComparisonConfig.yml"

with open(config["primary"]["file"], "r") as handle:
    config["primary"]["run"] = yaml.safe_load(handle)

with open(config["secondary"]["file"], "r") as handle:
    config["secondary"]["run"] = yaml.safe_load(handle)

if config["tertiary"]["file"]:
    with open(config["tertiary"]["file"],"r") as handle:
        config["tertiary"]["run"] = yaml.safe_load(handle)


module other_workflow:
    snakefile: "snakefile.smk"
    config: config["primary"]["run"]

use rule * from other_workflow as first_*

module another_workflow:
    snakefile: "snakefile.smk"
    config: config["secondary"]["run"]

use rule * from another_workflow as second_ *

if config["tertiary"]["file"]:

    module yetanother_worklow:
        snakefile: "snakefile.smk"
        config: config["tertiary"]["run"]

    use rule * from yetanother_worklow as third_ *


rule compareRuns:
    input:
        fileRun1 = rules.first_extractDESeqResult.output.result_table,
        fileRun2 = rules.second_extractDESeqResult.output.result_table,
        fileRun3 = branch(config["tertiary"]["file"], rules.third_extractDESeqResult.output.result_table, None)
    output:
        html = os.path.join(config["RUN_DIR"], "Comparison/DESeqResult_c{condition}_vs_b{baseline}.html"),
        svg = os.path.join(config["RUN_DIR"], "Comparison/DESeqResult_c{condition}_vs_b{baseline}.svg")
    run:
        from pyfunctions.plotting import compare_results
        to_compare = []
        names = [config[p]["name"] for p in ["primary", "secondary", "tertiary"]]
        for name, file in zip(names, input):
            to_compare.append((name, file))

        fig = compare_results(
            to_compare,sep="\t",condition=wildcards.condition,
            base=wildcards.baseline,padj_cutoff=0.05,lfc_cutoff=0.8,y=config["yName"]
        )
        fig.write_html(output.html)
        fig.write_image(output.svg)


rule all:
    default_target: True
    input:
        file = expand(rules.compareRuns.output.html, zip, condition=config["primary"]["run"]["conditions"], baseline=config["primary"]["run"]["baselines"])


