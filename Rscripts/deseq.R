suppressPackageStartupMessages({
  require(DESeq2)
  require(pheatmap)
})

countFile <- snakemake@input[["counts"]]
annotationFile <- snakemake@input[["annotation"]]
n_counts_output <- snakemake@output[["normalized_counts"]]
size_factors_output <- snakemake@output[["size_factors"]]

design_string <- snakemake@params[["design"]]

if (is.null(design_string) || design_string == "~" || design_string == "") {
    stop("Invalid or missing design formula passed from Snakemake. Expected something like '~ condition'")
}

formula <- as.formula(design_string)

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
norm_factors <- sizeFactors(dds)
write.table(norm_factors, size_factors_output, sep="\t", quote=FALSE)

write.table(normalized_counts, file=n_counts_output, sep="\t", quote=FALSE)

rld <- rlog(dds, blind=TRUE)
terms <- all.vars(formula)

# Take just the last term — typically the grouping variable
intgroup <- tail(terms, 1)
ntop <- min(500, nrow(rld))

pca_data <- plotPCA(rld, intgroup=terms, returnData=TRUE, ntop=500)
write.table(pca_data, snakemake@output[["pca_data"]], sep="\t")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
write.table(rld_cor, snakemake@output[["corr_data"]], sep="\t")
pheatmap(rld_cor, filename=snakemake@output[["heatmap"]])
png(snakemake@output[["dispersion_estimates"]])
plotDispEsts(dds)
dev.off()
save(dds, file=snakemake@output[["deseq_result"]])





