#!/usr/bin/env Rscript

################################################################
###### FUNCTIONS

ks.test.from.list <- function(geneid){
  
  dist_frag_gene <- strsplit(as.character(replicate[geneid,1]),",")
  dist_frag_gene <- as.numeric(unlist(dist_frag_gene))
  
  d1 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="two.sided")
  d2 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="less")
  d3 <- ks.test(dist_frag_gene,dist_frag_replicate, alternative="greater")
  
  out_two.sided <- paste(d1$statistic,d1$p.value,sep = ";")
  out_less <- paste(d2$statistic,d2$p.value,sep = ";")
  out_greater <- paste(d3$statistic,d3$p.value,sep = ";")
  output_line <- paste(gene, out_less, out_two.sided, out_greater, sep = " ")

}

##################################################################################################
######### import files
# try it: read.table(gzfile("/tmp/foo.csv.gz"))

## list of genes and their mapped read lengths
replicate<- read.table(args[1], header=F);
## list of genes
genes <- read.table(args[2], header = F)$V1
  
dist_frag_replicate <- strsplit(as.character(replicate[,1]),",")
dist_frag_replicate <- as.numeric(unlist(dist_frag_replicate))
  
for(i in 1:length(genes)){
  ks.test.from.list(as.character(genes[i]))
}




