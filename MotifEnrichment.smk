import yaml
import os
from Bio import SeqIO

configfile: "motifEnrichment.yml"

with open(config["DESeqConfig"]["file"], "r") as handle:
    config["DESeqConfig"]["run"] = yaml.safe_load(handle)


module other_workflow:
    snakefile: "snakefile.smk" if not config["multilevel"] else "multilevelDESeq.smk"
    config: config["DESeqConfig"]["run"]

use rule * from other_workflow

deseq_config = config["DESeqConfig"]["run"] if not config["multilevel"] else config["DESeqConfig"]["run"]["secondary"]


rule getFasta:
    input:
        gff = deseq_config["gff"],
        fasta = config["Fasta"]
    conda:
        "envs/bedtools.yml"
    output:
        fasta = temporary(os.path.join(config["RUN_DIR"],"fasta/sequences.fa")),
    shell:
        "bedtools getfasta -s -fi {input.fasta} -bed {input.gff} > {output.fasta}"

rule extractFasta:
    input:
        file = rules.extractDESeqResult.output.result_table if not config["multilevel"] else rules.second_extractDESeqResult.output.result_table,
        gff = deseq_config["gff"],
        fasta = rules.getFasta.output.fasta
    output:
        up_fasta = os.path.join(config["RUN_DIR"],"fasta/{condition}_vs_{baseline}_up.fa"),
        bg_fasta = os.path.join(config["RUN_DIR"],"fasta/{condition}_vs_{baseline}_bg.fa")
    run:
        import pandas as pd
        df = pd.read_csv(input.file, sep="\t")
        up = df[df["log2FoldChange"] >= config["motifEnrichment"]["lfc_cutoff"]]
        up = up[up["padj"] >= config["motifEnrichment"]["padj_cutoff"]]
        gff = pd.read_csv(input.gff, sep="\t",  names=["chr", "feature", "feat", "start", "end", "dot", "strand", "score", "add"], skiprows=1)
        gff["locus_tag"] = gff["add"].str.split("locus_tag=").str[-1]
        gff["up"] = gff["locus_tag"].isin(up.index)
        gff["match"] = gff["chr"] + ":" + (gff["start"] - 1).astype(str) + "-" + gff["end"].astype(str) + "(" + gff["strand"] + ")"
        up = gff[gff["up"]]
        bg = gff[~gff["up"]]
        ups = []
        bgs = []
        d= up["match"].tolist()
        for seq_record in SeqIO.parse(input.fasta, "fasta"):
            if str(seq_record.description) in d:
                ups.append(seq_record)
            else:
                bgs.append(seq_record)
        with open(output.up_fasta, "w") as handle:
            SeqIO.write(ups, handle, "fasta")

        with open(output.bg_fasta, "w") as handle:
            SeqIO.write(bgs, handle, "fasta")

rule memeMotif:
    input:
        up_fasta = rules.extractFasta.output.up_fasta,
        bg_fasta = rules.extractFasta.output.bg_fasta,
    conda:
        "envs/meme.yml"
    output:
        dir = directory(os.path.join(config["RUN_DIR"],"MEME/{condition}_vs_{baseline}_up/"))
    shell:
        "streme --p {input.up_fasta} --o {output.dir} --n {input.bg_fasta} --evalue"


rule allMotif:
    default_target: True
    input:
        file = expand(rules.memeMotif.output, zip, condition=deseq_config["conditions"], baseline=deseq_config["baselines"])
        #file = rules.getFasta.output
