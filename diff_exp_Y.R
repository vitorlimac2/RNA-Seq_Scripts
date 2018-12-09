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
my_quant <- "group_count_edited.txt"

featureCounts_file <- read.table(my_quant, header = F)

rownames(featureCounts_file) <- do.call(paste,c(featureCounts_file[c(1,2)],sep="  "))

cts <- featureCounts_file[,3:8]

colnames(cts) <- replicates
cts <- round(cts)
coldata <- data.frame(replicate=replicates, condition=c("Epi","Epi","Epi","Trypo","Trypo","Trypo"))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

dds <- dds[ rowSums(counts(dds)) > 1, ]

### Main result - DEG

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Trypo", "Epi"))

write.table(res, file="results_edited.csv", sep = "\t", quote = F)


#### Mean normalized count by condition

write.table(file="norm_counts_edited.txt",counts(dds, normalized = TRUE), sep = "\t", quote = F)


## 

tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj))
head(tab)
par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~padj))
lfc = 1.5
pval = 0.05
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("padj =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

