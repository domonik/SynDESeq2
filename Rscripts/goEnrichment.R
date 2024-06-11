suppressPackageStartupMessages({
  library(clusterProfiler)
})
package <- list.files(snakemake@input[["annotation_db"]])[1]
print(snakemake@input[["annotation_db"]])
library(basename(package), character.only = TRUE, lib.loc=snakemake@config[["Rlib"]])

log2cutoff <- snakemake@config[["log2FCCutOff"]]
padjcutoff <- snakemake@config[["pAdjCutOff"]]

minGSSize <- snakemake@config[["minGSSize"]]
maxGSSize <- snakemake@config[["maxGSSize"]]

defile <- snakemake@input[["deseq_results"]]
detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
has_X_column <- any(colnames(detable) == "X")
if (has_X_column) {
    # Remove unnamed column first
    # Set rownames from original first column header
    rownames(detable) <- detable[, 1]
    # Remove first row, now containing data
    detable <- detable[, -1 ]
}

up <- detable$log2FoldChange >= log2cutoff & detable$padj < padjcutoff
up <- detable[up, ]
down <- detable$log2FoldChange <= -log2cutoff & detable$padj < padjcutoff
down <- detable[down, ]
universe <- as.character(rownames(detable))
egoBPup <- enrichGO(gene = as.character(rownames(up)),
                    universe = universe,
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

ont <- "ALL"
ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))
op <- clusterProfiler:::get_GO_data(basename(package), ont, "GID")
TERMID2EXTID <- getFromNamespace("TERMID2EXTID", "DOSE")

filter_genes <- function(genes, valid_genes) {
  genes[genes %in% valid_genes]
}

qExtID2TermID <- TERMID2EXTID(as.character(summary$ID), op)
qExtID2TermID <- lapply(qExtID2TermID, filter_genes, valid_genes = universe)

summary$universeGeneID <- sapply(qExtID2TermID, function(x) paste(x, collapse = "/"))

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

qExtID2TermID <- TERMID2EXTID(as.character(summary$ID), op)
qExtID2TermID <- lapply(qExtID2TermID, filter_genes, valid_genes = universe)

summary$universeGeneID <- sapply(qExtID2TermID, function(x) paste(x, collapse = "/"))

write.table(summary , file = snakemake@output[["down"]], row.names=FALSE, sep="\t")

