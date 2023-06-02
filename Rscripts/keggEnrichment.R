suppressPackageStartupMessages({
  library(clusterProfiler, lib.loc = snakemake@config[["Rlib"]])
})
defile <- snakemake@input[["defile"]]



detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
detable <- detable[!grepl("UTR", rownames(detable)), ]


padjcutoff <- snakemake@config[["pAdjCutOff"]]
log2fccutoff <- snakemake@config[["log2FCCutOff"]]
org <- snakemake@config[["keggOrgID"]]

up <- detable$log2FoldChange >= log2fccutoff & detable$padj < padjcutoff
up <- detable[up, ]
down <- detable$log2FoldChange <= -log2fccutoff & detable$padj < padjcutoff
down <- detable[down, ]
ekeggup <- enrichKEGG(gene = as.character(rownames(up)),
                      universe=as.character(rownames(detable)),
                      organism=org,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
)
summary <- data.frame(ekeggup )
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 9, nrow = 0))
  x <- c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary, file = snakemake@output[["up"]], row.names = FALSE, sep = "\t")

ekeggdown <- enrichKEGG(gene = as.character(rownames(down)),
                        universe=as.character(rownames(detable)),
                        organism = org,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
summary <- data.frame(ekeggdown)
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 9, nrow = 0))
  x <- c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary, file = snakemake@output[["down"]], row.names=FALSE, sep="\t")