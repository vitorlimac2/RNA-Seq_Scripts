## Manipulate counts with deseq2

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library("DESeq2")
library("ggplot2")

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
                "ProR2",
                "HiSeq1",
                "HiSeq2")

plot.gene.ranking <- function(replicate, reference, norm_counts, rank = NULL, raw_counts = NULL, reference2 = NULL){
  
  rank = NULL
  if(is.null(rank)){
    target_genes <- rownames(norm_counts[order(norm_counts[,replicate], decreasing = T),])
  }else{
    target_genes <- rownames(norm_counts[order(norm_counts[,replicate], decreasing = T),])[1:rank]
    
  }
  plot(main = replicate,log10(norm_counts[target_genes,replicate]), 
       col="blue", 
       type="l",
       ylab="Counts (log10)",
       xlab="Gene Ranking", 
      log="x",
       ylim=c(1,7))
  if(!is.null(reference2)){
    points(log10(norm_counts_filtered[target_genes,reference2]), col="darkgray", pch ="*")
  }
  
  points(log10(norm_counts[target_genes,reference]), col="darkgreen", pch = "*")
  
  if(!is.null(raw_counts)){
    lines(log10(raw_counts_filtered[target_genes,replicate]), col="red")
  }
  
  if(is.null(reference2) & is.null(raw_counts)){
    legend("bottomleft",
           legend=c(reference, replicate),
           fill=c("darkgreen", "blue"),
           bty="n")
  }else if(!is.null(reference2) & !is.null(raw_counts)){
    legend("bottomleft",
           legend=c(reference, reference2, paste(replicate," (I)", sep=""),
                    paste(replicate," (II)",sep="")),
           fill=c("darkgreen","gray", "blue", "red"),
           bty="n")
  }else if(is.null(reference2) & !is.null(raw_counts)){
    legend("bottomleft",
           legend=c(reference, paste(replicate," (I)", sep=""),
                    paste(replicate, " (II)",sep="")),
           fill=c("darkgreen", "blue", "red"),
           bty="n")
  }else if(!is.null(reference2) & is.null(raw_counts)){
    legend("bottomleft",
           legend=c(reference, reference2, paste(replicate," (I)", sep="")),
           fill=c("darkgreen","gray", "blue"),
           bty="n")
  }
}

workdir <- "/home/vitor/Downloads"

setwd(workdir)

raw_counts <- "Counts_NoNormalized_Frac.txt"
norm_counts <- "Counts_Normalized_Frac.txt"

raw_counts <- read.table(raw_counts, header=FALSE, skip = 1,row.names = 1)
norm_counts <- read.table(norm_counts, header=FALSE, skip = 1, row.names = 1)

head(raw_counts)
head(norm_counts)
#raw_counts <- raw_counts[,3:(2+length(replicates))]
#norm_counts <- norm_counts[,3:(2+length(replicates))]

colnames(raw_counts) <- replicates
colnames(norm_counts) <- replicates

raw_counts_filtered <- raw_counts[rowMeans(raw_counts)>=1 & rowMeans(norm_counts) >= 1,]
norm_counts_filtered <- norm_counts[rowMeans(raw_counts)>=1 & rowMeans(norm_counts) >= 1,] 

replicate <- "ProR1"

reference <- "HiSeq1"
reference2 <- "HiSeq2"

target_genes <- rownames(norm_counts[order(norm_counts[,replicate], decreasing = T),])


plot.gene.ranking(replicate = replicate,
                  reference = reference,
                  norm_counts = norm_counts_filtered,
                  reference2 = reference2,
                  raw_counts = raw_counts_filtered)

identify(c(1:length(log10(norm_counts[target_genes,reference]))), 
         log10(norm_counts[target_genes,reference]), 
        target_genes, tolerance = 0.38, cex = 0.6, offset = 0.5)

identify(c(1:length(log10(norm_counts_filtered[target_genes,reference2]))),
         log10(norm_counts_filtered[target_genes,reference2]), 
         target_genes, tolerance = 0.38, cex = 0.6, offset = 0.5)
