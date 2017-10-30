#!/usr/bin/env Rscript

###########################################################################################
### Functions
ks.table <- function(tab_row){
  
  gene <- tab_row[1]
  dist_frag <- strsplit(tab_row[2],",")
  dist_frag <- as.numeric(unlist(dist_frag))
  
  d1 <- ks.test(rep_frags,dist_frag, alternative="two.sided")
  d2 <- ks.test(rep_frags,dist_frag, alternative="less")
  d3 <- ks.test(rep_frags,dist_frag, alternative="greater")
  
  out_two.sided <- paste(d1$statistic,d1$p.value,sep = ";")
  out_less <- paste(d2$statistic,d2$p.value,sep = ";")
  out_greater <- paste(d3$statistic,d3$p.value,sep = ";")
  
  output_line <- paste(gene, out_less, out_two.sided, out_greater, sep = " ")
  print(output_line,stdout())
  
}

###########################################################################################
###########################################################################################
######## ARGS

args <- commandArgs(trailingOnly=TRUE)


if(length(args)!=3){
  
  message("USAGE:\nRscript --vanilla KS_test_Dist_Frag.R repFrag mappedFrag numGenes\nOPTIONS:\n\trepFrag: File with two columns. 1st column contains read id. 2nd colunm contains read length.\n\tmappedFrag: File with two columns. 1st column contains gene id. 2nd column contains comma-separated list of fragment lengths.\n\tnumGenes: number of detected genes (number of lines of mappedFrag.)", call.=FALSE)
  stop("Missing options.")
}

replicate_read_lengths_file <- args[1]
inputFile <- args[2]
total_genes <- as.numeric(args[3])

rep_frags <- read.table(replicate_read_lengths_file, header = F)$V2

inputFile <- "/home/vitor/Proj_ProC_R/mappings_star/ProC1.Gene.DistFrag.DistMaps.ProteinCoding"
mapped_frags <- read.table(inputFile, header = F)[,1:2]

#######################################################################################################


#########################################################################################

## read each line

# while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
#     myVector <- strsplit(oneLine, " ")
#     gene <- myVector[[1]][1]
#     dist_frag <- strsplit(myVector[[1]][2],",")
#     
#     dist_frag <- as.numeric(unlist(dist_frag))
#     
#     d1 <- ks.test(rep_frags,dist_frag, alternative="two.sided")
#     d2 <- ks.test(rep_frags,dist_frag, alternative="less")
#     d3 <- ks.test(rep_frags,dist_frag, alternative="greater")
#     
#     out_two.sided <- paste(d1$statistic,d1$p.value,sep = ";")
#     out_less <- paste(d2$statistic,d2$p.value,sep = ";")
#     out_greater <- paste(d3$statistic,d3$p.value,sep = ";")
#     
#     output_line <- paste(gene, out_less, out_two.sided, out_greater, sep = " ")
#     
#     myFrac <- i*100/total_genes; 
#     
#     i = i + 1;
#     
#     my_log <- paste("Running ", myFrac,"%", sep="")
#     
#     print(my_log, stderr())
#     print(output_line,stdout())
#     
# }
# 
# close(con)