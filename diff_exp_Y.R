## quantification with DESeq2

source("https://bioconductor.org/biocLite.R")
#library(bioconductor)
biocLite("DESeq2")

a
library(DESeq2)
library(ggplot2)

replicates <- c("Epi1",
                "Epi2",
                "Epi3",
                "Trypo1",
                "Trypo2",
                "Trypo3")

## calculate normalized factors

#my_quant <- "allReplicates_output_protein.coding"
my_quant <- "group_count.txt"

setwd("/home/vitor/Projects/Paralog_Quantification_Y/")

featureCounts_file <- read.table(my_quant, header = F)

rownames(featureCounts_file) <- do.call(paste,c(featureCounts_file[c(1,2)],sep="  "))

head(featureCounts_file)
  
cts <- featureCounts_file[,3:8]

head(cts,2)
colnames(cts) <- replicates

head(cts,1)

cts <- round(cts)
coldata <- data.frame(replicate=replicates, condition=c("Epi","Epi","Epi","Trypo","Trypo","Trypo"))

coldata

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

head(cts)

nrow(counts(dds))

dds <- dds[ rowSums(counts(dds)) > 1, ]

nrow(counts(dds))


######################### The variance stabilizing transformation and the rlog
## plot the standard deviation of each row (genes) against the mean
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

# And for logarithm-transformed counts:
  
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)




vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)


library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)


sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = "condition")

pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()


mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed()


mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed()


dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast=c("Epi", "Trypo"))

res

write.csv(res, file="results.csv")

###############################################
###############################################
###############################################


dat <- stack(as.data.frame(log10(counts(dds))))

ggplot(dat) + geom_boxplot(aes(x = ind, y = values)) + labs(x="Replicates", y="Raw counts (log10)")

cpms <- counts(dds)

head(cpms)

cpms[,1] <- counts(dds)[,1]*1000000/sum(counts(dds)[,1])
cpms[,2] <- counts(dds)[,2]*1000000/sum(counts(dds)[,2])
cpms[,3] <- counts(dds)[,3]*1000000/sum(counts(dds)[,3])
cpms[,4] <- counts(dds)[,4]*1000000/sum(counts(dds)[,4])
cpms[,5] <- counts(dds)[,5]*1000000/sum(counts(dds)[,5])
cpms[,6] <- counts(dds)[,6]*1000000/sum(counts(dds)[,6])
cpms[,7] <- counts(dds)[,7]*1000000/sum(counts(dds)[,7])
cpms[,8] <- counts(dds)[,8]*1000000/sum(counts(dds)[,8])
cpms[,9] <- counts(dds)[,9]*1000000/sum(counts(dds)[,9])
cpms[,10] <- counts(dds)[,10]*1000000/sum(counts(dds)[,10])
cpms[,11] <- counts(dds)[,11]*1000000/sum(counts(dds)[,11])
  
head(cpms)

dat <- stack(as.data.frame(log10(cpms)))

ggplot(dat) + geom_boxplot(aes(x = ind, y = values)) + labs(x="Replicates", y="Counts per million (log10)")



cts
dds <- estimateSizeFactors(dds)
rld <-rlog(dds)

## rlog transformation

plotPCA(rld, intgroup=c("condition","type")) + 
  geom_text(aes(label = replicates), position = position_nudge(y = 1)) + 
  ggtitle("Read Counts - rlog transformation") +
  theme(plot.title = element_text(hjust=0.5,size=14))

### LOG2 transformation + normalized by sizeFactors

se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))

plotPCA( DESeqTransform(se), intgroup=c("condition","type")) + 
  geom_text(aes(label = replicates), position = position_nudge(y = 1)) + 
  ggtitle("Read Counts - normalized by sizeFactors (log2)") +
  theme(plot.title = element_text(hjust=0.5,size=14))


## normalized by sizeFactors

se <- SummarizedExperiment(counts(dds, normalized=TRUE),
                           colData=colData(dds))

plotPCA( DESeqTransform(se), intgroup=c("condition", "type")) + 
  geom_text(aes(label = replicates), position = position_nudge(y = 1)) + 
  ggtitle("Read Counts - normalized by sizeFactors (w/o log2)") +
  theme(plot.title = element_text(hjust=0.5,size=14))

## raw

se <- SummarizedExperiment(counts(dds),
                           colData=colData(dds))
plotPCA( DESeqTransform(se), intgroup=c("condition", "type")) + 
  geom_text(aes(label = replicates), position = position_nudge(y = 1)) + 
  ggtitle("Read Counts - without transformations") +
  theme(plot.title = element_text(hjust=0.5,size=14))

dds <- DESeq(dds)

res <- results(dds)
head(res[order(res$log2FoldChange),])
