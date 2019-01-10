#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
inputFile <- args[1]
cts <- read.table(inputFile, row.names = 1, header = T)
replicates <- c("ProC1", "ProC2", "ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2","Hiseq1","Hiseq2")
cts <- cts[rowSums(cts) > 1,]
cts <- round(cts)
coldata <- data.frame(replicate=replicates, 
                      condition=c("ProC","ProC","ProC","ProC","ProC","ProC","ProC","ProC","ProC",
                                                        "ProR", "ProR",
                                                        "HiSeq","HiSeq"),
                      initRNA = c(2,2,2,2,2,.2,.2,.2,.2,2,2,2,2))
library("DESeq2",quietly = T, warn.conflicts = F)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition+initRNA)
dds <- estimateSizeFactors(dds)
cts <- as.data.frame(counts(dds, normalized = TRUE))
write.table(cts, paste(inputFile,"SizeFactorEstimated", sep = "."), sep = "\t", quote = F, col.names = F)

