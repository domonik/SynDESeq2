import os
import re

import pandas as pd
import requests
from requests.adapters import Retry, HTTPAdapter

rule copy_config_file:
    output:
        file = os.path.join(config["RUN_DIR"], "config.yml")
    run:
        import yaml
        print(config)
        with open(output.file, "w") as handle:
            yaml.dump(config, handle)


def annotation_from_structured_counts(file, sep, column_mapping):
    df = pd.read_csv(file, sep=sep, index_col=0)
    data = [[] for _ in range(len(column_mapping))]
    names = []
    for col in df.columns:
        names.append(col)
        name = col.split("_")
        for idx, sublist in enumerate(data):
            sublist.append(name[idx])
    annotation_dict = {
        key: data[idx] for idx, key in enumerate(column_mapping)
    }
    annotation = pd.DataFrame(annotation_dict)
    annotation.index = names
    return annotation


rule generateAnnotationFromCounts:
    input:
        counts = config["CountFile"],
        cfg = rules.copy_config_file.output.file
    output:
        annotation = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/Annotation.csv")
    run:
        annotation = annotation_from_structured_counts(
            input.counts, sep=config["initial_sep"], column_mapping=config["column_mapping"]
        )
        annotation.to_csv(output.annotation, sep="\t")


def drop_unused_columns(counts, annotation, drop_conditions, sep, blacklist):

    counts = pd.read_csv(counts, sep=sep, index_col=0)
    counts.index = counts.index.str.replace("'", "")
    annotation = pd.read_csv(annotation, sep="\t", index_col=0)
    for key, item_list in drop_conditions.items():
        for item in item_list:
            counts = counts.loc[:, annotation[key] != item]
            annotation = annotation[annotation[key] != item]
    if blacklist:
        for item in blacklist:
            counts = counts.loc[:, annotation.index != item]
            annotation = annotation[annotation.index != item]

    return counts, annotation


rule dropUnusedSamples:
    input:
        counts = config["CountFile"],
        annotation = rules.generateAnnotationFromCounts.output.annotation,
    output:
        counts = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DroppedCounts.tsv"),
        annotation = os.path.join(config["RUN_DIR"], "PipelineData/IntermediateData/DroppedAnnotation.tsv")
    run:
        counts, annotation = drop_unused_columns(
            input.counts, input.annotation, config["drop_condition"], sep=config["initial_sep"], blacklist=config["blacklisted_samples"]
        )

        counts.to_csv(output.counts, sep="\t")
        annotation.to_csv(output.annotation, sep="\t")


ORGANISMID = f"{config['genus']}.{config['species']}_taxid{config['tax_id']}"


def download_organism_go_terms(tax_id, outfile):
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    #url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cgo&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cxref_ensembl%2Cgo_id%2Cft_transmem%2Cft_intramem&format=tsv&query=%28%28taxonomy_id%3A{tax_id}%29%29&size=500"
    #purl = "https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cxref_ensembl&format=tsv&query=%28%28taxonomy_id%3A9606%29%29&size=500"
    progress = 0
    print("Starting to Download GO Terms")
    with open(outfile, "w") as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'Downloaded {progress} / {total}')
    if progress == 0:
        raise ValueError("The outfile is empty.")
    return 0


rule downloadOrganismGOTerms:
    output:
        go_terms = os.path.join(config["GOTermDir"],  ORGANISMID + ".tsv")
    run:
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


