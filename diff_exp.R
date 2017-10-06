## quantification with DESeq2

#source("https://bioconductor.org/biocLite.R")
#library(bioconductor)
#biocLite("DESeq2")

library(DESeq2)
library(ggplot2)
## read complete table

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

## calculate normalized factors

#my_quant <- "allReplicates_output_protein.coding"
my_quant <- "allReplicates_output_allGenes"

featureCounts_file <- read.table(my_quant, row.names = 1, header = T)

head(featureCounts_file)

cts <- featureCounts_file[,6:16]

colnames(cts) <- replicates

head(cts)
coldata <- read.table("colData", row.names = 1, header = T)

coldata

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

head(cts)



dds

nrow(counts(dds))

my_limit <- 0


dds <- dds[counts(dds)[,1] > my_limit & counts(dds)[,2] > my_limit & counts(dds)[,3] > my_limit &
             counts(dds)[,4] > my_limit & counts(dds)[,5] > my_limit & counts(dds)[,6] > my_limit &
             counts(dds)[,7] > my_limit & counts(dds)[,8] > my_limit & counts(dds)[,9] > my_limit &
           counts(dds)[,10] > my_limit & counts(dds)[,11] > my_limit, ]
nrow(counts(dds))

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
