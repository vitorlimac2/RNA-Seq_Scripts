#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
################################################################
###### FUNCTIONS

ks.test.from.list <- function(geneid){
  
  dist_frag_gene <- strsplit(as.character(levels(replicate.GeneDistFrag[geneid,1]))[replicate.GeneDistFrag[geneid,1]],",")
  dist_frag_gene <- as.numeric(unlist(dist_frag_gene))

  d1 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="two.sided")
  d2 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="less")
  d3 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="greater")
  
  out_two.sided <- paste(d1$statistic,d1$p.value,sep = ";")
  out_less <- paste(d2$statistic,d2$p.value,sep = ";")
  out_greater <- paste(d3$statistic,d3$p.value,sep = ";")
  output_line <- paste(geneid, out_less, out_two.sided, out_greater, sep = " ")
  
  return(output_line)

}

##################################################################################################
######### import files
# try it: read.table(gzfile("/tmp/foo.csv.gz"))

## list of genes and their mapped read lengths
replicate.GeneDistFrag <- read.table(args[1], header=F, row.names=1);

#replicate<- read.table("/media/vitor/Seagate Expansion Drive/Thesis/DEG_ILB/ILB_9577.GeneDistFrag", header=F, row.names = 1);

## list of genes
genes <- read.table(args[2], header = F)$V1

dist_frag_replicate <- strsplit(as.character(levels(replicate.GeneDistFrag[,1])[replicate.GeneDistFrag[,1]]),",")
dist_frag_replicate <- as.numeric(unlist(dist_frag_replicate))

genes <- as.character(levels(genes))[genes]
  
for(i in 1:length(genes)){
  print(ks.test.from.list(genes[i]), quote=F)
}
