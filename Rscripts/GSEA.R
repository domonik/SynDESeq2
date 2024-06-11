suppressPackageStartupMessages({
  library(clusterProfiler)
})
package <- list.files(snakemake@input[["annotation_db"]])[1]
library(basename(package), character.only = TRUE, lib.loc=snakemake@config[["Rlib"]])

table <- snakemake@input[["deseq_results"]]
table <- read.table(table, sep="\t", header=TRUE)
table <- table[!is.na(table$log2FoldChange), ]
df_sorted <- table[order(table$log2FoldChange, decreasing = TRUE), ]
geneList <- df_sorted$log2FoldChange
names(geneList) <- rownames(df_sorted)
ego <- gseGO(geneList = geneList,
              OrgDb = basename(package),
              ont = "ALL",
              keyType = "GID",
              minGSSize = 30,
              maxGSSize = 500,
              pvalueCutoff = 0.4,
              verbose = FALSE,
)
summary <- data.frame(ego)
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 12, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
  colnames(df) <- x
  summary <- df
}

write.table(summary , file = snakemake@output[["enriched"]], row.names=FALSE, sep="\t")

