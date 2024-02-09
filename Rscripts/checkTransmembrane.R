suppressPackageStartupMessages({
  library(clusterProfiler, lib.loc = snakemake@config[["Rlib"]])
})

defile <- snakemake@input[["deseq_results"]]
term2name <- snakemake@input[["term2name"]]
log2cutoff <- snakemake@config[["log2FCCutOff"]]
padjcutoff <- snakemake@config[["pAdjCutOff"]]


term2name <- read.table(term2name, sep="\t", header=TRUE)

detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
up <- detable$log2FoldChange >= log2cutoff & detable$padj < padjcutoff
up <- detable[up, ]

res <- enricher(
  as.character(rownames(up)),
  universe=rownames(detable),
  TERM2GENE=term2name,
  minGSSize = 3,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  qvalueCutoff = 0.9,
)
summary <- data.frame(res )
write.table(summary , file = snakemake@output[["up"]], row.names=FALSE, sep="\t")

