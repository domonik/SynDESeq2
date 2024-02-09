suppressPackageStartupMessages({
  require(DESeq2)
  require(pheatmap)
})

countFile <- snakemake@input[["counts"]]
annotationFile <- snakemake@input[["annotation"]]
n_counts_output <- snakemake@output[["normalized_counts"]]

formula <- as.formula(snakemake@params[["design"]])

print(formula)

countData <- as.matrix(read.table(countFile, header = T,  sep="\t", row.names=1, check.names = F, comment.char = ""))
annotationData <- as.matrix(read.table(annotationFile, header = T, row.names = 1, sep="\t", check.names = F, comment.char = ""))

if (snakemake@params[["use_spike_ins"]]){
  spike_ins <- grepl(snakemake@params[["spike_in_pattern"]], rownames(countData))
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
  dds <- estimateSizeFactors(dds, controlGenes=spike_ins)

  dds <- dds[!spike_ins,]

} else if (snakemake@params[["use_housekeeping"]]) {
  housekeeping <- read.table(snakemake@input[["housekeeping"]], header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "")
  spike_ins <- row.names(countData) %in% row.names(housekeeping)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
  dds <- estimateSizeFactors(dds, controlGenes=spike_ins)
} else {
  countData <- as.matrix(countData)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
}
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file=n_counts_output, sep="\t", quote=FALSE)

rld <- rlog(dds, blind=TRUE)
pca_data <- plotPCA(rld, intgroup="Fraction", returnData=TRUE, ntop=500)
write.csv(pca_data, snakemake@output[["pca_data"]], sep="\t")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, filename=snakemake@output[["heatmap"]])
png(snakemake@output[["dispersion_estimates"]])
plotDispEsts(dds)
dev.off()
save(dds, file=snakemake@output[["deseq_result"]])





