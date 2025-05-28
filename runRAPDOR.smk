import yaml


configfile: "runRAPDOR.yml"

with open(config["DESeqConfig"]["file"], "r") as handle:
    config["DESeqConfig"]["run"] = yaml.safe_load(handle)

module other_workflow:
    snakefile: "snakefile.smk"
    config: config["DESeqConfig"]["run"]

use rule * from other_workflow as first_*

rule prepareForRAPDOR:
    input:
        normalized_counts = rules.first_runDESeq.output.normalized_counts,
        design = rules.first_dropUnusedSamples.output.annotation,
        annotation_file = branch(config["Annotation"]["file"], config["Annotation"]["file"], None)
    output:
        rapdor_design = os.path.join(config["RUN_DIR"],"initData/init_design.tsv"),
        rapdor_counts = os.path.join(config["RUN_DIR"],"initData/normalized_counts.tsv")
    run:
        import pandas as pd
        from pyfunctions.run_rapdor import generateRAPDORInput
        df = pd.read_csv(input.normalized_counts, sep="\t", index_col=0)
        design = pd.read_csv(input.design, sep="\t", index_col=0)
        design = design.rename(
            {
                config["treatment"]["col"]: "Treatment",
                config["fraction"]["col"]: "Fraction",
                config["replicate"]["col"]: "Replicate"
            },
            axis=1
        )
        design.index.name = "Name"
        if config["Annotation"]["file"]:
            annotation = pd.read_csv(input.annotation_file, sep=config["Annotation"]["sep"], index_col=0)
            df = pd.merge(df,annotation, left_index=True, right_index=True, how="left")
        if "gene_name" in df:
            df = df.rename({"gene_name": "Gene"}, axis=1)
        df.to_csv(output.rapdor_counts, sep="\t")
        design.to_csv(output.rapdor_design, sep="\t")


rule runRAPDOR:
    input:
        counts = rules.prepareForRAPDOR.output.rapdor_counts,
        design = rules.prepareForRAPDOR.output.rapdor_design
    output:
        json = os.path.join(config["RUN_DIR"],"json/RAPDOR.json"),
    run:
        import pandas as pd
        from RAPDOR.datastructures import RAPDORData
        df = pd.read_csv(input.counts, sep="\t")
        design = pd.read_csv(input.design, sep="\t")

        rapdordata = RAPDORData(df=df,design=design, control=config["treatment"]["control"], measure="Counts", measure_type="Gene",
            min_replicates=3)
        rapdordata.normalize_and_get_distances(method="Jensen-Shannon-Distance")
        rapdordata.calc_all_scores()
        rapdordata.rank_table(["ANOSIM R", "Mean Distance"],ascending=(False, False))
        print(rapdordata.df["Rank"])
        rapdordata.to_json(output.json)


rule all:
    default_target: True
    input:
        file = rules.runRAPDOR.output

