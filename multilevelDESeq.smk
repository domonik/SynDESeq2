include: "rules/enrich.smk"

configfile: "multilevelConfig.yml"

rule all:
    input:
        file =   os.path.join(config["RUN_DIR"], "PipelineData/Multilevel/houskeeping_genes2.tsv")

rule runOnLeastDiffRegGenes:
    input:
        files = expand(rules.extractDESeqResult.output.result_table,  zip, condition=config["conditions"], baseline=config["baselines"])
    output:
        file =   os.path.join(config["RUN_DIR"], "PipelineData/Multilevel/intersectingGenes.tsv")
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

def extract_intersect(file, sep: str = "\t"):
    import pandas as pd
    df = pd.read_csv(file, sep=sep)
    return df.index.tolist()

rule multilevelDESeq:
    input:
        files=expand(rules.extractDESeqResult.output.result_table,zip,condition=config["conditions"],baseline=config[
            "baselines"])

    params:
        use_housekeeping =  config["use-housekeeping"],
        use_spike_ins = config["use-spike-ins"],
        design = config["design"],
        housekeeping = extract_intersect(rules.runOnLeastDiffRegGenes.output.file)
    output:
        file = os.path.join(config["RUN_DIR"], "PipelineData/Multilevel/houskeeping_genes.tsv")
    run:
        print("Success")