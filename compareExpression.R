### Ler arquivo

wd <- "/media/vitor/Seagate Expansion Drive/Thesis/experiment.CompareNormalization/"
setwd(wd)

inputFile <- "ProC1.vs.ProR1.Samples"

cts <- read.table(inputFile, row.names = 1)[,1:4]

cts <- cts[rowSums(cts) > 10,]

cts <- round(cts)

replicates <- c("ProC1_raw", "ProR1_raw", 
                "ProC1_norm", "ProR1_norm")

coldata <- data.frame(replicate=replicates, 
                      condition=c("ProC","ProR","ProC","ProR"),
                      normalization = c("raw","raw", "norm", "norm"))

library("DESeq2",quietly = T, warn.conflicts = F)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition+normalization)
                      
dds <- DESeq(dds)




vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, fitType = "local")

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$normalization, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition","normalization")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



plotPCA(vsd, intgroup=c("condition", "normalization"), )                      
vsd

vds
## Plotar

head(f)
library(ggplot2)

f$Read_Lengths <- factor(f$Read_Lengths, levels = c("Control", "ShortReads", "LongReads"), labels = c("Control", "ShortReads", "LongReads"))

cts <- round(f)

ggplot(f,aes(x=log(ProR1_raw), y = log(ProC1_raw), color=Read_Lengths)) + geom_point() + geom_smooth()

ggplot(f,aes(x=log(ProR1_norm), y = log(ProC1_norm), color=Read_Lengths)) + geom_point()+ geom_smooth()


cor.test(f$ProC1_norm, f$ProR1_norm, method=c("pearson", "kendall", "spearman"))

library("ggpubr")
