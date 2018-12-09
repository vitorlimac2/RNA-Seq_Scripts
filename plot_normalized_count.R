wd <- "/media/vitor/Seagate Expansion Drive/Thesis/"

setwd(wd)

library("DESeq2")
library("ggplot2")

replicates <- c("ProR1",
                "ProR2",
                "HiSeq1",
                "HiSeq2")
my_quant <- "ProR1_ProR2_HiSeq1_HiSeq2_RawCounts.Frac.Multi"

featureCounts_file <- read.table(my_quant, header = F, row.names = 1)

cts <- featureCounts_file

colnames(cts) <- replicates
nrow(cts)
cts <- cts[cts$ProR1 > 0 & cts$ProR2 > 0 & cts$HiSeq1 > 0 & cts$HiSeq2 > 0,]
nrow(cts)
cts <- round(cts)


coldata <- data.frame(replicate=replicates, condition=c("Ion","Ion","Illumina","Illumina"))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
#dds <- DESeq(dds)
#res <- results(dds, contrast=c("condition","Illumina", "Ion"))
cts_norm <- as.data.frame(counts(dds, normalized = TRUE))


## 

plot(log10(cts_norm[order(cts_norm$ProR1, decreasing = T),]$HiSeq1),
     col="burlywood4",
     xlab="Ranking", 
     ylab="Normalized count (log10)",
     log="x",
     ylim= c(0,7),
     xlim=c(1000,10000),
     lwd=4,
     pch="*")

lines(col="black",
      log10(cts_norm[order(cts_norm$ProR1, decreasing = T),]$ProR1), 
      lwd=3)      


identify(c(1:length(log10(cts_norm[order(cts_norm$ProR1, decreasing = T),]$HiSeq1))),
         log10(cts_norm[order(cts_norm$ProR1, decreasing = T),]$HiSeq1),
         rownames(cts_norm[order(cts_norm$ProR1, decreasing = T),]), cex=0.8, col="red")

points(x = c(258,300),
       y = c(log10(cts_norm[order(cts_norm$ProR1, decreasing = T),][258,"HiSeq1"]),log10(cts_norm[order(cts_norm$ProR1, decreasing = T),][300,"HiSeq1"])), col="red", pch=5)

legend("bottomleft", 
       legend=c("Ion Torrent Proton","Illumina HiSeq"), 
       col=c("black", "burlywood4"), inset = 0.05,box.lty=0,
       cex = 1, 
       lty = 1,
       lwd=3)
