#!/usr/bin/env Rscript

###########################################################################################
### Functions
ks.table <- function(tab_row){
  
  gene <- tab_row[1]

  #dist_frag <- strsplit(as.character(tab_row[2]),",")
  dist_frag <- strsplit(as.character(tab_row[2]),",")
  dist_frag <- as.numeric(unlist(dist_frag))
  
  d1 <- ks.test(dist_frag,rep_frags, alternative="two.sided")
  d2 <- ks.test(dist_frag,rep_frags, alternative="less")
  d3 <- ks.test(dist_frag,rep_frags, alternative="greater")
  
  out_two.sided <- paste(d1$statistic,d1$p.value,sep = ";")
  out_less <- paste(d2$statistic,d2$p.value,sep = ";")
  out_greater <- paste(d3$statistic,d3$p.value,sep = ";")
  output_line <- paste(gene, out_less, out_two.sided, out_greater, sep = " ")
  
  return(output_line)
  
}

###########################################################################################
###########################################################################################
######## ARGS

args <- commandArgs(trailingOnly=TRUE)


if(length(args)!=3){
  message("USAGE:\nRscript --vanilla KS_test_Dist_Frag.R repFrag mappedFrag numCores\nOPTIONS:\n\trepFrag: File with two columns. 1st column contains read id. 2nd colunm contains read length.\n\tmappedFrag: File with two columns. 1st column contains gene id. 2nd column contains comma-separated list of fragment lengths.\n\tnumCores: number of threads", call.=FALSE)
  stop("Missing options.")
}

# uncomment this at the end
replicate_read_lengths_file <- args[1]
inputFile <- args[2]
numCores <- as.numeric(args[3])

## comment this at the end
# setwd("/media/vitor/Seagate Expansion Drive/Thesis/DEG_ILB/")
# replicate_read_lengths_file <- "Test.MappedReadInfo"
# inputFile <- "Test.GeneDistFrag"
# numCores <- as.numeric("2")



#### ONLY FOR TEST ##################
#replicate_read_lengths_file <- "/home/vitor/Proj_ProC_R/reads/trimmed/ProC_1.fastq.trimmed.fq.read_lengths.txt"
#inputFile <- "/home/vitor/Proj_ProC_R/mappings_star/ProC1.Gene.DistFrag.DistMaps.ProteinCoding.test"
#total_genes <- 5
#numCores <- 2
######################################
rep_frags <- read.table(replicate_read_lengths_file, header = F)
rep_frags <- rep_frags[rep_frags$V3 == 1,]$V2
mapped_frags <- read.table(inputFile, header = F)
total_genes <- nrow(mapped_frags)
#######################################################################################################
### MULTI-THREAD

library(parallel,quietly = T, verbose = F)
## Create cluster
clus <- makeCluster(numCores, type="FORK")

## The clusterExport() function exports an object to each node, 
## enabling them to work parallely. The use of it, as it can be 
# appreciated, is extremely simple: you need to pass the variable 
# name/s in a character vector (or a single string, as in this case)
clusterExport(clus,"rep_frags")

aa <- parApply(clus,mapped_frags,1, ks.table)
output_file <- paste(inputFile,".ks.output",sep="")
write(aa,output_file)
stopCluster(clus)
gc()

#########################################################################################

# ## TEST
# 
# # "ENSMUSG00000050737 0.409950145409223;0.00012788302329305 0.409950145409223;0.000255766046585504 0.0447652679684254;0.89862359276185"
# 
# gene <- "ENSMUSG00000050737"
# 
# #dist_frag <- strsplit(as.character(tab_row[2]),",")
# dist_frag <- strsplit(as.character(mapped_frags[mapped_frags$V1 == "ENSMUSG00000050737",2]),",")
# dist_frag <- as.numeric(unlist(dist_frag))
# 
# ks.test(dist_frag,rep_frags, alternative="two.sided")
# ks.test(dist_frag,rep_frags, alternative="less")
# ks.test(dist_frag,rep_frags, alternative="greater")
