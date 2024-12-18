#Load libraries
library(edgeR)
#library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)
library(DESeq2)

#Before running the code, change working directory to the directory of the
#files being imported

#PART 1

#Import RNA count and sample table
counts <- read.table(file = "E-MTAB-2523.counts.txt",
                          header = TRUE,
                          as.is = TRUE,
                          sep = "\t",
                          row.names = 1)

samples <- read.table(file = "E-MTAB-2523_sample table.txt",
                          header = TRUE,
                          as.is = TRUE,
                          sep = "\t",
                          row.names = 1)

#Data filtering, based on mean log2 cmp
meanLog2CPM <- rowMeans(log2(cpm(counts) + 1))
counts <- counts[meanLog2CPM >1, ]

#PART 2

#Extract data from sample table, this data is later used for design of dds
#object
individual <- samples$individual
disease <- samples$disease

#Create dds object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  design = ~ individual + disease,
  colData = data.frame(
    individual = individual,
    disease = disease))

dds <- DESeq(dds)
res <- results(dds)

#FDR and log2 fold change
Filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]





