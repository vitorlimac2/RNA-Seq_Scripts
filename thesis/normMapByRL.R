#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=1){
  message("USAGE:\nRscript --vanilla normMapByRL.R <file>.MappedReadInfo.\nOPTIONS:\n\t<file>.MappedReadInfo: 4-column output file of Obtain_mapped_read_information.sh. 1st = read_id; 2nd = read length; 3rd = number of mappings; 4st = mappings status (1= mapped; 0 = unmapped)", call.=FALSE)
  stop("Missing options.")
}

## Test
## input file
## READ | LENGTH | STATUS | GENE

## Read input file

## File with read, it length and mapping status

inputFile <- args[1]

f1 <- read.table(inputFile, fill = T, sep="\t");

read_col <- 1;
length_col <- 2;
n_mappings <- 3;
mapping_status <- 4;

## Calculate bin vector of assigned read lengths
length_bins <- pretty(f1[f1$V4 == 1,]$V2);

length_bins[1] <- min(f1$V2);
length_bins[length(length_bins)] <- max(f1$V2);

mapping_rate_vector <- vector(mode="numeric",length=length(length_bins));
mapped_read_vector <- vector(mode="numeric",length=length(length_bins));
unmapped_read_vector<- vector(mode="numeric",length=length(length_bins));

## Calculate mapping rate for each bin 

for(i in seq_len(length(length_bins)-1)){
  x <- i;
  y <- i+1;
  
  if(y != length(length_bins)){
    mapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y] & f1$V4==1,]);
    unmapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 < length_bins[y] & f1$V4==0,]);
    }else{
    mapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y] & f1$V4==1,]);
    unmapped_read_vector[y] <- nrow(f1[f1$V2 >= length_bins[x] & f1$V2 <= length_bins[y] & f1$V4==0,]);
  }
}
## overall mapping rate
mapping_rate_vector <- mapped_read_vector/(mapped_read_vector+unmapped_read_vector)*100;

## Calculate median mapping rate
m <- median(mapping_rate_vector,na.rm = T);

## Calculate normalization factor for each bin
normalization_factor <- m/mapping_rate_vector;

## Assign the normalization factor to each read

temp_file <- f1[f1$V4!=0,];

rm("f1");

temp_file$V5 <- 0

for(i in seq_len(length(length_bins)-1)){
  x <- i;
  y <- i+1;
  
  if(y != length(length_bins)){
    temp_file[temp_file$V2 >= length_bins[x] & temp_file$V2 < length_bins[y],]$V5 <- normalization_factor[y];
  }else{
    temp_file[temp_file$V2 >= length_bins[x] & temp_file$V2 <= length_bins[y],]$V5 <- normalization_factor[y];
  }
}

## Print the new gene count with normalization factor

write.table(temp_file,paste(inputFile,"ReadCountNormalizedByRL",sep = "."),sep = "\t", quote = F, row.names = F, col.names = F);

rm("temp_file")
