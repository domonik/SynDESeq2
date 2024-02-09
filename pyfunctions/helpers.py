import pandas as pd
import requests
import re
import requests
from requests.adapters import HTTPAdapter, Retry
from sklearn.cluster import DBSCAN
import numpy as np


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
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cft_transmem%2Cft_intramem&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
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
    return 0


def cluster_go_enrich(files, enrich_file, config):
    enrich_df = pd.read_csv(enrich_file, sep="\t")
    if len(enrich_df) == 0:
        enrich_df["Cluster"] = None
        return enrich_df

    rows = []
    clusters = []
    for file in files:
        df = pd.read_csv(file, sep="\t", index_col=0)
        if np.isnan(df.iloc[0, 0]):
            continue
        distance_matrix = 1 - np.asarray(df.iloc[:, :])
        clustering = DBSCAN(metric="precomputed", eps=config["cluster-eps"], min_samples=2)
        fitted = clustering.fit_predict(distance_matrix)
        rows += list(df.columns)
        clusters += list(fitted)

    clusters = pd.DataFrame(
        {
            "ID": rows,
            "Cluster": clusters
        }

    )
    df = pd.merge(enrich_df, clusters, on="ID")
    return df


def reduce_cluster(df: pd.DataFrame, method: str = "strlen", ascending: bool = True):
    non_clustered = df[df.Cluster == -1]
    df = df[~(df.Cluster == -1)]
    if method == "strlen":
        df["sort_column"] = df["Description"].str.len()
        ascending = ascending
    elif method in df.columns:
        df["sort_column"] = df[method]
        ascending = ascending
    else:
        raise NotImplementedError(f"Method '{method}' not implemented")
    df = df.sort_values("sort_column", ascending=ascending)
    df = df.drop_duplicates(["ONTOLOGY", "Cluster"], keep='first')
    df = pd.concat((df, non_clustered), axis=0)
    df = df.drop("sort_column", axis=1)
    return df

if __name__ == '__main__':
    file = "../LightPuromycinMinus/PipelineData/GOEnrichment/SemanticSimilarity_ontCC_up_cM_vs_bTC.tsv"
    cluster_go_enrich(file)

