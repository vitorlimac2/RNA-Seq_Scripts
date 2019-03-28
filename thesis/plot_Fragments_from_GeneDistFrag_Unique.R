#!/usr/bin/env Rscript

################################################################
###### FUNCTIONS


fragment.plot <- function(geneid){
  
  dist_frag_gene <- strsplit(as.character(replicate[geneid,1]),",")
  dist_frag_gene <- as.numeric(unlist(dist_frag_gene))
  
  my_breaks <- pretty(range(c(dist_frag_gene, dist_frag_replicate)), n=5)
  
  h_1 <- hist(dist_frag_replicate, breaks=my_breaks, plot=FALSE)
  h_2 <- hist(dist_frag_gene, breaks=my_breaks, plot=FALSE)
  
  max_y <- max(h_1$counts*100/length(dist_frag_replicate), h_2$counts*100/length(dist_frag_gene))
  
  #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  id <- paste(replicate_name,geneid)
  l_info <- paste("L(",test.ks.qvalue[geneid,]$less,
                     ",",
                     test.ks.qvalue[geneid,]$less_distance,
                     ")",
                     sep="")
  t_info <- paste("\nT(",test.ks.qvalue[geneid,]$twosided,
                     ",",
                     test.ks.qvalue[geneid,]$twosided_distance,
                     ")",
                     sep="")
  g_info <- paste("G(",test.ks.qvalue[geneid,]$greater,
                     ",",
                     test.ks.qvalue[geneid,]$greater_distance,
                     ")",
                     sep="")
  info <- paste(l_info,t_info,g_info, sep=";")
  
  my_title <- paste(id,info)
  
  png(paste(replicate_name,geneid,"fragmentation","png", sep="."))
  
  plot(main = my_title, 
       h_1$mids, 
       h_1$counts*100/length(dist_frag_replicate), 
       type="l", 
       ylab = "%", 
       xlab = "Read length (bp)", 
       col="blue", 
       ylim=c(0,max_y), 
       cex.main=0.8)
  
  lines(h_2$mids, h_2$counts*100/length(dist_frag_gene), col="red")
  dev.off()
}

##################################################################################################
######### import files
# try it: read.table(gzfile("/tmp/foo.csv.gz"))
setwd("/media/vitor/Seagate Expansion Drive/Thesis/DEG")

replicates <- c("ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2")

for(j in 1:length(replicates)){
  replicate_name <- replicates[j]
  replicate<- read.table(paste(replicate_name, ".GeneDistFrag_Unique",sep=""), row.names = 1)
  
  test.ks.qvalue <- read.table(paste(replicate_name, ".top200.ExpressionDiff.KS.Qvalues", sep =""), row.names = 1)
  
  colnames(test.ks.qvalue) <- c("less_distance", "less", "twosided_distance", "twosided", "greater_distance", "greater")
  
  genes <- read.table(paste(replicate_name, ".top200.ExpressionDiff.KS.Qvalues", sep = ""))$V1
  
  dist_frag_replicate <- strsplit(as.character(replicate[,1]),",")
  dist_frag_replicate <- as.numeric(unlist(dist_frag_replicate))
  
  for(i in 1:length(genes)){
    fragment.plot(as.character(genes[i]))
  }
}



