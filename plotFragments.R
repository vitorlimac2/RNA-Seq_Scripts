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
  
  plot(main = geneid, h_1$mids, h_1$counts*100/sum(length(replicateFrags)), type="l", ylab = "%", xlab = "Read length (bp)", col="blue", ylim=c(0,max_y), xlim=c(35,250)) 
  
  lines(h_2$mids, h_2$counts*100/sum(length(dist_frag_gene)), col="red")
}

##################################################################################################
######### import files
# try it: read.table(gzfile("/tmp/foo.csv.gz"))
replicateFrags <- read.table("ProC1.ReadLengths")$V2
mappings_table <- read.table("ProC1.Gene.DistFrag.DistMaps.ProteinCoding", row.names = 1)
#pvalues_table <- read.table("ProC1.Gene.DistFrag.DistMaps.ProteinCoding.ks.output.Pvalues", header=F,row.names = 1)

file_input <- "ProC1.Gene.DistFrag.DistMaps.ProteinCoding.ks.output.Pvalues.Qvalues"
qvalues_table <- read.table(file_input, header = T, row.names = 1)
detach(qvalues_table)
attach(qvalues_table)

x <- rownames(qvalues_table[less <= 0.01 & two_sided <= 0.01 & greater > 0.01,])


output_file <- paste(file_input,".Less", sep = "")
write(x,output_file, sep = "\n")

x <- rownames(qvalues_table[less > 0.01 & two_sided <= 0.01 & greater > 0.01,])
output_file <- paste(file_input,".2sided", sep = "")
write(x,output_file,sep="\n")

x <- rownames(qvalues_table[less > 0.01 & two_sided <= 0.01 & greater <= 0.01,])
output_file <- paste(file_input,".Greater", sep = "")
write(x,output_file,sep="\n")

x <- rownames(qvalues_table[less <= 0.01 & two_sided <= 0.01 & greater <= 0.01,])
output_file <- paste(file_input,".All", sep = "")
write(x,output_file,sep="\n")




qvalues_table <- qvalues_table[order(qvalues_table$less),]
qvalues_matrix <- data.matrix(qvalues_table[less <= 0.05 | two_sided <= 0.05 | greater <= 0.05,])

qheatmap <- heatmap(qvalues_matrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))

heatmap(qvalues_matrix,Colv=NA,scale='none')
