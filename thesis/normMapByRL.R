#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=1){
  message("USAGE:\nRscript --vanilla normMapByRL.R <file>.MappedReadInfo.\nOPTIONS:\n\t<file>.MappedReadInfo: 4-column output file of Obtain_mapped_read_information.sh. 1st = read_id; 2nd = read length; 3rd = number of mappings; 4st = mappings status (1= mapped; 0 = unmapped)", call.=FALSE)
  stop("Missing options.")
}

## input file
## READ | LENGTH | NUMBER_OF_MAPPINGS | STATUS 

inputFile <- args[1]

f1 <- read.table(inputFile, sep="\t")

nrow(f1)
length_bins <- pretty(f1$V2, n=5)

#length_bins <- seq(min(f1$V2),max(f1$V2),by = 20)

length_bins[1] <- min(f1$V2)
length_bins[length(length_bins)] <- max(f1$V2);

length_bins

mapping_rate_vector <- vector(mode="numeric",length=length(length_bins));
mapped_read_vector <- vector(mode="numeric",length=length(length_bins));
unmapped_read_vector<- vector(mode="numeric",length=length(length_bins));

## Calculate mapping rate for each bin 

for(i in seq_len(length(length_bins)-1)){
  x <- i;
  y <- i+1;
  
  if(y != length(length_bins)){
    mapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y] & f1$V4==1,]);
    unmapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y] & f1$V4!=1,]);
  }else{
    mapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y] & f1$V4==1,]);
    unmapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y] & f1$V4!=1,]);
  }
}

## overall mapping rate
mapping_rate_vector <- mapped_read_vector*100/(mapped_read_vector+unmapped_read_vector);
mapping_rate_vector
x <- mapping_rate_vector[mapping_rate_vector!=Inf]

## Calculate median mapping rate
m <- median(x,na.rm = T);

m

## Calculate normalization factor for each bin
normalization_factor <- m/mapping_rate_vector;
normalization_factor
## Assign the normalization factor to each read

f1$V5 <- 0

for(i in seq_len(length(length_bins)-1)){
  x <- i;
  y <- i+1;
  
  if(y != length(length_bins)){
    if(normalization_factor[y]=="Inf" | normalization_factor[y]=="NaN"){
      f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y],]$V5 <- 0;
    }else{
      f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y],]$V5 <- normalization_factor[y];
    }
    
  }else if(normalization_factor[y]=="Inf" | normalization_factor[y]=="NaN"){
    f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y],]$V5 <- 0;
  }else{
    f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y],]$V5 <- normalization_factor[y];
    
  }
}

summary(f1$V5)

## Print the new gene count with normalization factor

write.table(f1,paste(inputFile,"ReadCountNormalizedByRL",sep = "."),sep = "\t", quote = F, row.names = F, col.names = F);

