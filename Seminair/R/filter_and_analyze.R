#' Filter RNA-seq Counts and Identify DEGs
#'
#' This function filters lowly expressed genes from RNA-seq count data using EdgeR's `filterByExpr` method
#' and identifies differentially expressed genes (DEGs) using DESeq2.
#' Filtered counts and DEGs are saved to the specified output directory.
#'
#' @param counts_file Path to the RNA-seq counts file.
#' @param samples_file Path to the sample metadata file.
#' @param output_dir Directory to save results. Defaults to "output".
#' @return A list containing filtered counts and significant DEGs.
#' @export
filter_and_analyze <- function(counts_file, samples_file, output_dir = "output") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  library(edgeR)
  library(DESeq2)
  counts <- read.table(file = counts_file, header = TRUE, sep = "\t", row.names = 1)
  samples <- read.table(file = samples_file, header = TRUE, sep = "\t", row.names = 1)

  group <- factor(samples$disease)
  DEG <- DGEList(counts = counts, group = group)
  filter <- filterByExpr(DEG)
  DEG <- DEG[filter, , keep.lib.sizes = FALSE]

  write.table(DEG$counts, file = file.path(output_dir, "filtered_counts.txt"), sep = "\t", quote = FALSE)

  dds <- DESeqDataSetFromMatrix(DEG$counts, colData = data.frame(disease = group), design = ~ disease)
  dds <- DESeq(dds)
  res <- results(dds)

  significant_degs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
  write.table(as.data.frame(significant_degs), file = file.path(output_dir, "significant_degs.txt"), sep = "\t", quote = FALSE)

  list(filtered_counts = DEG$counts, significant_degs = significant_degs)
}

#usage
counts_file <- "input/E-MTAB-2523.counts.txt"
samples_file <- "input/E-MTAB-2523_sample table.txt"
output_dir <- "output"
results <- filter_and_analyze(counts_file, samples_file, output_dir)
