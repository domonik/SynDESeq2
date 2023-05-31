suppressPackageStartupMessages({
  require(DESeq2)
  require(pheatmap)
})

countFile <- snakemake@input[["counts"]]
annotationFile <- snakemake@input[["annotation"]]
formula <- as.formula(snakemake@config[["design"]])

print(formula)

countData <- as.matrix(read.table(countFile, header = T,  sep="\t", row.names=1, check.names = F, comment.char = ""))
annotationData <- as.matrix(read.table(annotationFile, header = T, row.names = 1, sep="\t", check.names = F, comment.char = ""))

if (snakemake@config[["use-spike-ins"]]){
  spike_ins <- grepl(snakemake@config[["spike-in-pattern"]], rownames(countData))
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
  dds <- estimateSizeFactors(dds, controlGenes=spike_ins)

  dds <- dds[!spike_ins,]

} else if (snakemake@config[["use-housekeeping"]]) {
  housekeeping <- read.table(snakemake@config[["housekeeping"]], header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "")
  spike_ins <- row.names(countData) %in% row.names(housekeeping)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
  dds <- estimateSizeFactors(dds, controlGenes=spike_ins)
} else {
  countData <- as.matrix(countData)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)
}
dds <- DESeq(dds)

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





