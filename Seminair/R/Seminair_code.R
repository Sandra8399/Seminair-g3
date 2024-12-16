#Load libraries
library(edgeR)
#library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)
library(DESeq2)

#Before running the code, change working directory to the directory of the
#files being imported

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



