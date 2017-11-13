#!/usr/bin/env Rscript

################################################################
###### FUNCTIONS


fragment.plot <- function(geneid, geneMappings, replicateFrags){
  
  # geneId = string
  # geneMappings = file *.Gene.DistFrag.DistMaps.ProteinCoding
  # replicateFrags = file *.ReadLengths
  
  dist_frag_gene <- strsplit(as.character(geneMappings[geneid,1]),",")
  dist_frag_gene <- as.numeric(unlist(dist_frag_gene))
  
  my_breaks <- pretty(range(c(dist_frag_gene, replicateFrags)), n=10)
  
  h_1 <- hist(replicateFrags, breaks=my_breaks, plot=FALSE)
  h_2 <- hist(dist_frag_gene, breaks=my_breaks, plot=FALSE)
  
  max_y <- max(h_1$counts*100/sum(length(replicateFrags)), h_2$counts*100/sum(length(dist_frag_gene)))
  
  #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  plot(h_1$mids, h_1$counts*100/sum(length(replicateFrags)), type="l", ylab = "%", xlab = "Read length (bp)", col="blue", ylim=c(0,max_y), xlim=c(35,250)) 
  
  lines(h_2$mids, h_2$counts*100/sum(length(dist_frag_gene)), col="red")
}

##################################################################################################
######### import files
replicateFrags <- read.table("ProC1.ReadLengths")$V2
mappings_table <- read.table("ProC1.Gene.DistFrag.DistMaps.ProteinCoding", row.names = 1)
pvalues_table <- read.table("ProC1.Gene.DistFrag.DistMaps.ProteinCoding.ks.output.Pvalues", header=F,row.names = 1)
qvalues_table <- read.table("ProC1.Gene.DistFrag.DistMaps.ProteinCoding.ks.output.Pvalues.Qvalues", header = T, row.names = 1)


nrow(qvalues_table[less <= 0.05 & two_sided > 0.05 & greater > 0.05,])
nrow(qvalues_table[two_sided <= 0.05,])
nrow(qvalues_table[greater <= 0.05,])
nrow(qvalues_table[less <= 0.05 & two_sided <= 0.05 & greater > 0.05,])
nrow(qvalues_table[greater <= 0.05 & two_sided <= 0.05 & less > 0.05,])
nrow(qvalues_table[greater <= 0.05 & two_sided <= 0.05 & less <= 0.05,])

nrow(qvalues_table[less <= 0.05,])*100/nrow(qvalues_table)
nrow(qvalues_table[two_sided <= 0.05,])*100/nrow(qvalues_table)
nrow(qvalues_table[greater <= 0.05,])*100/nrow(qvalues_table)
nrow(qvalues_table[less <= 0.05 & two_sided <= 0.05 & greater > 0.05,])*100/nrow(qvalues_table)
nrow(qvalues_table[greater <= 0.05 & two_sided <= 0.05 & less > 0.05,])*100/nrow(qvalues_table)
nrow(qvalues_table[greater <= 0.05 & two_sided <= 0.05 & less <= 0.05,])*100/nrow(qvalues_table)
