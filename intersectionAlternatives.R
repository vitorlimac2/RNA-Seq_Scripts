###
### Replicate names

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
                "ProR2")

## alternatives

input_file_terminal <- ".Gene.DistFrag.DistMaps.ProteinCoding.ks.output.Pvalues.Qvalues.Greater"

for(i in seq(1,length(replicates)-1)){
  for(j in seq(i+1,length(replicates))){
  
  rep1 <- read.table(paste(replicates[i],input_file_terminal,sep=""))$V1
  rep2 <- read.table(paste(replicates[j],input_file_terminal,sep=""))$V1
  
  intsec <- length(intersect(rep1,rep2))
  un <- length(union(rep1,rep2))
  
  print(paste(replicates[i],replicates[j],un))
  }
}
