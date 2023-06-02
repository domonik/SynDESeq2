suppressPackageStartupMessages({
  library(clusterProfiler, lib.loc = snakemake@config[["Rlib"]])
})
package <- list.files(snakemake@input[["annotation_db"]])[1]
library(basename(package), character.only = TRUE)

log2cutoff <- snakemake@config[["log2FCCutOff"]]
padjcutoff <- snakemake@config[["pAdjCutOff"]]

minGSSize <- snakemake@config[["minGSSize"]]
maxGSSize <- snakemake@config[["maxGSSize"]]

defile <- snakemake@input[["deseq_results"]]
detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)

up <- detable$log2FoldChange >= log2cutoff & detable$padj < padjcutoff
up <- detable[up, ]
down <- detable$log2FoldChange <= -log2cutoff & detable$padj < padjcutoff
down <- detable[down, ]

egoBPup <- enrichGO(gene = as.character(rownames(up)),
                    universe = as.character(rownames(detable)),
                    OrgDb = basename(package),
                    keyType = "GID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = FALSE)
summary <- data.frame(egoBPup )
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary , file = snakemake@output[["up"]], row.names=FALSE, sep="\t")

egoBPdown <- enrichGO(gene = as.character(rownames(down)),
                    universe = as.character(rownames(detable)),
                    OrgDb = basename(package),
                    keyType = "GID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = TRUE)

summary <- data.frame(egoBPdown )
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary , file = snakemake@output[["down"]], row.names=FALSE, sep="\t")

