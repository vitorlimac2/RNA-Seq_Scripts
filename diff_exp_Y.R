## quantification with DESeq2

#source("https://bioconductor.org/biocLite.R")
#library(bioconductor)
#biocLite("DESeq2")

library("DESeq2")
library("ggplot2")

replicates <- c("Epi1",
                "Epi2",
                "Epi3",
                "Trypo1",
                "Trypo2",
                "Trypo3")

setwd("/media/vitor/Seagate Expansion Drive/Luciana/")
my_quant <- "membrane_counts_group_filtered_edited.tsv"

featureCounts_file <- read.table(my_quant, header = F)

rownames(featureCounts_file) <- do.call(paste,c(featureCounts_file[c(1,2)],sep="  "))

cts <- featureCounts_file[,3:8]

colnames(cts) <- replicates
cts <- round(cts)
coldata <- data.frame(replicate=replicates, condition=c("Epi","Epi","Epi","Trypo","Trypo","Trypo"))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

dds <- dds[ rowSums(counts(dds)) > 10, ]

### Main result - DEG

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Trypo", "Epi"), pAdjustMethod = "fdr")

write.table(res, file="results_filtered_edited.csv", sep = "\t", quote = F)


#### Mean normalized count by condition

write.table(file="norm_counts_filtered_edited.txt",counts(dds, normalized = TRUE), sep = "\t", quote = F)


##  volcano plot

tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj))
head(tab)
par(mar = c(5, 5, 3, 3))
plot(tab, pch = 16, cex = 0.6, xlab = "Fold-change (log2)", ylab = "Adjusted p-value (-log10)")
lfc = 1.5
pval = 0.05
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("padj =", pval), side = 4, at = -log10(pval), cex = 0.8, line = -4, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

### TRANSFORMING

ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

##### HEATMAP

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition")])

f <- read.table("median_counts_membrane_components.tsv", header = T, row.names = 1)

fn <- data.frame(epi=f$epi+0.01, trypo=f$trypo+0.01)
rownames(fn) <- rownames(f)

pheatmap(fn, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_row = f$component)


pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=coldata$replicates) 

### 

sampleDists <- dist(t(assay(vsd)[]))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
