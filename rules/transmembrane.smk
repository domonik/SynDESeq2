import os
import pandas as pd

include: "setup.smk"
include: "deseq.smk"


ORGANISMID = f"{config['genus']}.{config['species']}_taxid{config['tax_id']}"


rule prepareTransmembraneProteins:
    input:
        uniprot_table = rules.downloadOrganismGOTerms.output.go_terms,
    output:
        file = os.path.join(config["GOTermDir"],  ORGANISMID + "_TransmembraneTable.tsv")
    run:
        df = pd.read_csv(input.uniprot_table, sep="\t")
        df = df.rename({"Gene Names (ordered locus)": "locus_tag"},axis=1)
        df = df[~df["locus_tag"].isna()]
        df["Transmembrane"] = ~df["Transmembrane"].isna()
        df["Transmembrane"] = df["Transmembrane"].map({True: "Transmembrane", False: "No Transmembrane"})
        df = df[["Transmembrane", "locus_tag"]]
        df.to_csv(output.file, sep="\t", index=False)


rule transmembraneEnrichment:
    input:
        term2name = rules.prepareTransmembraneProteins.output.file,
        deseq_results = rules.extractDESeqResult.output.result_table
    conda:
        "../envs/REnvironment.yml"
    output:
        up = os.path.join(config["RUN_DIR"],"PipelineData/Enrichment/Transmembrane_up_c{condition}_vs_b{baseline}.tsv")
    script:
        "../Rscripts/checkTransmembrane.R"


rule downloadSynProteinDistribution:
    output:
        soluble = os.path.join(config["GOTermDir"], "SolubleProteins.xlsx") ,
        membrane = os.path.join(config["GOTermDir"], "MembraneProteins.xlsx") ,
    shell:
        """
        wget --no-check-certificate --keep-session-cookies --save-cookies=cookies --load-cookies=cookies --user-agent="Mozilla" -O {output.soluble} https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_003.xlsx
        wget --no-check-certificate --keep-session-cookies --save-cookies=cookies --load-cookies=cookies --user-agent="Mozilla" -O {output.membrane} https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_002.xlsx
        """


rule convertSynProteinDistribution:
    input:
        soluble = "initialData/SolubleProteins.xlsx",
        membrane = "initialData/MembraneProteins.xlsx"
    output:
        joined = os.path.join(config["RUN_DIR"],"PipelineData/Proteomics/JoinedProteomics.tsv")
    run:
        df = pd.read_excel(input.soluble, skiprows=2)
        df = df[["Ctrl1_1", "Ctrl1_2", "Ctrl1_3",]]
        print(df.iloc[2])


rule joinmRNALentgh:
    input:
        gff = config["gff"],
        deseq_result = rules.extractDESeqResult.output.result_table
    output:
        file = os.path.join(config["RUN_DIR"],"PipelineData/Lentgh/Correlation_c{condition}_vs_b{baseline}.tsv")
    run:
        from scipy.stats import spearmanr
        cols = ["chr", "feature_type", "feature_type2", "start", "end", "score", "strand", "score2", "add_info"]
        df = pd.read_csv(input.gff, sep="\t", comment="#", names=cols)
        df["locus_tag"] = df["add_info"].str.split("locus_tag=").str[-1].str.split(";").str[0]
        df["length"] = df["end"] -df["start"]
        df.index = df.locus_tag
        df = df[["length"]]

        deseq_res = pd.read_csv(input.deseq_result,sep="\t")

        deseq_res = deseq_res.join(df)
        statistic, pvalue = spearmanr(deseq_res["log2FoldChange"], deseq_res["length"], nan_policy="omit")
        data = pd.DataFrame(
            {
                "x": ["log2FoldChange"],
                "y": ["length"],
                "SpearmanR": [statistic],
                "p-Value": [pvalue]

            }
        )
        data.to_csv(output.file, sep="\t", index=False)
