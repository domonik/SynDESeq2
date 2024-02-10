import os

import pandas as pd

rule copy_config_file:
    output:
        file = os.path.join(config["RUN_DIR"], "config.yml")
    run:
        import yaml
        print(config)
        with open(output.file, "w") as handle:
            yaml.dump(config, handle)

rule generateAnnotationFromCounts:
    input:
        counts = config["CountFile"],
        cfg = rules.copy_config_file.output.file
    output:
        annotation = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/Annotation.csv")
    run:
        from pyfunctions.helpers import annotation_from_structured_counts
        annotation = annotation_from_structured_counts(
            input.counts, sep=config["initial_sep"], column_mapping=config["column_mapping"]
        )
        annotation.to_csv(output.annotation, sep="\t")


rule dropUnusedSamples:
    input:
        counts = config["CountFile"],
        annotation = rules.generateAnnotationFromCounts.output.annotation,
    output:
        counts = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DroppedCounts.tsv"),
        annotation = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DroppedAnnotation.tsv")
    run:
        from pyfunctions.helpers import drop_unused_columns
        counts, annotation = drop_unused_columns(
            input.counts, input.annotation, config["drop_condition"], sep=config["initial_sep"], blacklist=config["blacklisted_samples"]
        )

        counts.to_csv(output.counts, sep="\t")
        annotation.to_csv(output.annotation, sep="\t")


ORGANISMID = f"{config['genus']}.{config['species']}_taxid{config['tax_id']}"

rule downloadOrganismGOTerms:
    output:
        go_terms = os.path.join(config["GOTermDir"],  ORGANISMID + ".tsv")
    run:
        from pyfunctions.helpers import download_organism_go_terms
        download_organism_go_terms(config["tax_id"], output.go_terms)

rule prepareOrgGOTerms:
    input:
        uniprotgo = rules.downloadOrganismGOTerms.output.go_terms
    output:
        go_terms = os.path.join(config["GOTermDir"], "all_terms_" + ORGANISMID + ".tsv"),
        symbols = os.path.join(config["GOTermDir"], "symbols_" + ORGANISMID + ".tsv"),
    run:
        df = pd.read_csv(input.uniprotgo, sep="\t")
        df = df.rename({"Gene Names (ordered locus)": "locus_tag"}, axis=1)
        df = df[~df["locus_tag"].isna()]
        symbols = df[["Entry", "locus_tag"]]
        df["GOTerm"] = df["Gene Ontology IDs"].str.split('; | |/')
        df = df[["locus_tag", "GOTerm"]]
        df = df.explode("locus_tag").explode("GOTerm")
        df = df[~df["GOTerm"].isna()]
        df.to_csv(output.go_terms, sep="\t", index=False)
        symbols.to_csv(output.symbols, sep="\t", index=False)

rule generateOrgDB:
    input:
        symbols = rules.prepareOrgGOTerms.output.symbols,
        go_terms = rules.prepareOrgGOTerms.output.go_terms
    params:
        species=config["species"],
        genus=config["genus"],
    conda:
        "../envs/REnvironment.yml"
    output:
        annotation_db = directory(
            os.path.join(config["GOTermDir"], "AnnotationDBs", ORGANISMID)
        ),
        finished_file = temporary(
            os.path.join(config["GOTermDir"], "finished_" + ORGANISMID)
        )
    script:
        "../Rscripts/createAnnotation.R"


