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

my_quant <- "raw_counts_protein_coding.txt"

work_dir <- "/home/vitor/Proj_ProC_R/mappings_star/"

cts <- read.table(paste(work_dir,my_quant, sep=""), row.names = 1, header = F)

head(cts)

head(cts)

colnames(cts) <- replicates

head(cts)

conditions <- c("C","C","C","C","C","C","C","C","C","R","R")
initialRNA <- c("2","2","2","2","2","0.2","0.2","0.2","0.2","2","2")
insertSizes <- c("A","B","A","B","A","B","A","B","A","B","B")

#conditions <- c("C","C","C","C","C","C","C","C","C")
#initialRNA <- c("2","2","2","2","2","0.2","0.2","0.2","0.2")
#insertSizes <- c("A","B","A","B","A","B","A","B","A")


colData <- data.frame(condition = conditions, initialRNA = initialRNA, insertSize = insertSizes)

colData

rownames(colData) <- replicates

colData

keep <- rowMeans(cts) > 2

cts <- cts[keep,]

head(cts)



cts <- round(cts)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~initialRNA + condition)


dds <- estimateSizeFactors(dds)
rld <-rlog(dds)

## rlog transformation

plotPCA(rld, intgroup=c("condition","initialRNA", "insertSize")) + 
  geom_text(aes(label = replicates), position = position_nudge(y = 1)) + 
  ggtitle("Read Counts - rlog transformation") +
  theme(plot.title = element_text(hjust=0.5,size=14))


  cts

  #########################################################################################
dds <- DESeq(dds)

res <- results(dds)
head(res[order(res$log2FoldChange),])