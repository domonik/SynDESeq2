from DEplots.volcano import volcano_from_deseq_result
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input.deseq_result, sep="\t")
    add_data = pd.read_csv(snakemake.input.uniprot_table, sep="\t")
    add_data = add_data.set_index(snakemake.config["id_column"])
    df = df.join(add_data, how="left")
    highlight = {}
    if snakemake.config["highlight"]:
        for key, value in snakemake.config["highlight"].items():
            color, names = value
            print(names)
            mask = df.index.str.contains(names)
            to_highlight = df[mask]
            if len(to_highlight) == 0:
                to_highlight = df[df[snakemake.config["name_column"]].str.contains(names) == True]
                assert len(to_highlight) > 0, f"Nothing to highlight for {key}-{names}"
            highlight[key] = (color, to_highlight.index)

    fig = volcano_from_deseq_result(
        deseq_result=df,
        name_col=snakemake.config["name_column"],
        highlight=highlight,
        lfc_cut_off=snakemake.config["log2FCCutOff"],
        padj_cutoff=snakemake.config["pAdjCutOff"],
        condition_name=snakemake.wildcards["condition"],
        base_name=snakemake.wildcards["baseline"],

    )
    fig.write_html(snakemake.output.html)
    fig.write_image(snakemake.output.svg)
    fig.write_json(snakemake.output.json)
