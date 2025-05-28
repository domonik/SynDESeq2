suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(tidyr)
  library(purrr)
})
gseaScores <- getFromNamespace("gseaScores", "DOSE")

gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

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
              pvalueCutoff = 0.05,
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

m <- nrow(summary)
geneSetID <- 1:m

if (dim(summary)[1] == 0){
  x <- c("x", "runningScore", "position", "ymin", "ymax", "geneList", "gene", "Description")
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df) <- x
  gsdata <- df
} else {
  gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = ego))
}
write.table(gsdata , file = snakemake@output[["gsdata"]], row.names=FALSE, sep="\t")
