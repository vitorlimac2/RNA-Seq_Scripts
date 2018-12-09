#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly=TRUE)
setwd("/media/vitor/Seagate Expansion Drive/Thesis/")
#inputFile <- args[1]
inputFile <- "All.RawCount.Frac.Unique.txt"
cts <- read.table(inputFile, row.names = 1)
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
#dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "initRNA"))
res <- results(dds)
res
write.table(cts, paste(inputFile,"SizeFactorEstimated", sep = "."), sep = "\t", quote = F, col.names = F)
######################################################################
######################################################################

inputFile <- "All.NormalizedCount.Frac.Unique.txt"
cts <- read.table(inputFile, row.names = 1)
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
#dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "initRNA"), )

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$initRNA, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

####################################################################################
####################################################################################

inputFile <- "All.RawCount.Frac.Unique.txt"
cts <- read.table(inputFile, row.names = 1)
replicates <- c("ProC1", 
                "ProC2", 
                "ProC3",
                "ProC4",
                "ProC5",
                "ProC6",
                "ProC7",
                "ProC8",
                "ProC9",
                "ProR1",
                "ProR2")
cts <- cts[rowSums(cts) > 10,]
cts <- round(cts)
coldata <- data.frame(replicate=replicates, 
                      condition=c("Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "Chemical",
                                  "RnaseIII", 
                                  "RnaseIII"),
                      initRNA = c(2,2,2,2,2,.2,.2,.2,.2,2,2))
library("DESeq2")
head(cts)
cts <- cts[,1:11]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition + initRNA)
#dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "initRNA"))

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$initRNA, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

cts_n <- as.data.frame(counts(dds, ))