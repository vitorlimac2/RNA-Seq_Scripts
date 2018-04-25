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
  
  plot(main = geneid, h_1$mids, h_1$counts*100/sum(length(replicateFrags)), type="l", ylab = "%", xlab = "Read length (bp)", col="blue", ylim=c(0,max_y), xlim=c(120,210)) 
  
  lines(h_2$mids, h_2$counts*100/sum(length(dist_frag_gene)), col="red")
}

#####################################################
#### MAIN
#################

replicateFrags <- read.table("ProC2.ReadLengths")$V2
mappings_table <- read.table("ProC2.Gene.DistFrag.DistMaps.ProteinCoding", row.names = 1)

fragment.plot("ENSG00000000938",mappings_table, replicateFrags)