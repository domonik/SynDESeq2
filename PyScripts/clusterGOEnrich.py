import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

def cluster_go_enrich(files, enrich_file, config):
    enrich_df = pd.read_csv(enrich_file, sep="\t")
    if len(enrich_df) == 0:
        enrich_df["Cluster"] = None
        return enrich_df

    rows = []
    clusters = []
    for file in files:
        df = pd.read_csv(file, sep="\t", index_col=0)
        df_c =  df.dropna(axis=1, how="all").dropna(axis=0, how="all")
        if len(df_c) == 0:
            continue
        distance_matrix = 1 - np.asarray(df_c.iloc[:, :])
        clustering = DBSCAN(metric="precomputed", eps=config["cluster-eps"], min_samples=2)
        fitted = clustering.fit_predict(distance_matrix)
        all_entities = list(df.columns)
        kept_entities = list(df_c.columns)
        entity_to_cluster = dict(zip(kept_entities, fitted))

        # Build clusters list with O(1) lookup per entity
        clusters += [entity_to_cluster.get(entity, -1) for entity in all_entities]
        rows += all_entities

    clusters = pd.DataFrame(
        {
            "ID": rows,
            "Cluster": clusters
        }

    )
    df = pd.merge(enrich_df, clusters, on="ID")
    return df




if __name__ == '__main__':


    df = cluster_go_enrich(snakemake.input.semsim, snakemake.input.enrichment, config=snakemake.config)
    df.to_csv(snakemake.output.table, sep="\t", index=False)
