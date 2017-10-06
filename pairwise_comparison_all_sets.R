## 
setwd("~/PRJ.SRP064142/QUANTIFY_REPLICATES")
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

comparisons <- c("ProC1_ProC2",
                   "ProC1_ProC3",
                   "ProC1_ProC4",
                   "ProC1_ProC5",
                   "ProC1_ProC6",
                   "ProC1_ProC7",
                   "ProC1_ProC8",
                   "ProC1_ProC9",
                   "ProC2_ProC3",
                   "ProC2_ProC4",
                   "ProC2_ProC5",
                   "ProC2_ProC6",
                   "ProC2_ProC7",
                   "ProC2_ProC8",
                   "ProC2_ProC9",
                   "ProC3_ProC4",
                   "ProC3_ProC5",
                   "ProC3_ProC6",
                   "ProC3_ProC7",
                   "ProC3_ProC8",
                   "ProC3_ProC9",
                   "ProC4_ProC5",
                   "ProC4_ProC6",
                   "ProC4_ProC7",
                   "ProC4_ProC8",
                   "ProC4_ProC9",
                   "ProC5_ProC6",
                   "ProC5_ProC7",
                   "ProC5_ProC8",
                   "ProC5_ProC9",
                   "ProC6_ProC7",
                   "ProC6_ProC8",
                   "ProC6_ProC9",
                   "ProC7_ProC8",
                   "ProC7_ProC9",
                   "ProC8_ProC9",
                   "ProR1_ProC1",
                   "ProR1_ProC2",
                   "ProR1_ProC3",
                   "ProR1_ProC4",
                   "ProR1_ProC5",
                   "ProR1_ProC6",
                   "ProR1_ProC7",
                   "ProR1_ProC8",
                   "ProR1_ProC9",
                   "ProR1_ProR2",
                   "ProR2_ProC1",
                   "ProR2_ProC2",
                   "ProR2_ProC3",
                   "ProR2_ProC4",
                   "ProR2_ProC5",
                   "ProR2_ProC6",
                   "ProR2_ProC7",
                   "ProR2_ProC8",
                   "ProR2_ProC9")

num_comparisons <- length(comparisons)

file_terminal <- "_gene_median_cts_trim.nt_raw.nt"

comparison_list <- list()

for(i in replicates){
  name_comp <- paste(i,file_terminal,sep="")
  comparison_list[[i]] <- read.table(name_comp)
  
  print(i)
  print(head(comparison_list[[i]]))
  
}

##################################################################
################################# Plot read counts
limit_c <- 0
library(vioplot)
library(beanplot)

criteria_diff <- log10(3)

for(i in comparisons){
  x <- strsplit(i,"_");
  rep1 <- x[[1]][1]
  rep2 <- x[[1]][2]
  head(comparison_list[[rep1]])
  head(comparison_list[[rep2]])
  
  comparison <- cbind(comparison_list[[rep1]]$V2, comparison_list[[rep1]]$V3, comparison_list[[rep2]]$V3,
                      comparison_list[[rep1]]$V4, comparison_list[[rep2]]$V4)
  rownames(comparison) <- comparison_list[[rep1]]$V1
  head(comparison)
  
  mydata_count <- comparison
  
  mydata_count[,2] <- comparison[,2]*1000000/sum(comparison[,2])
  mydata_count[,3] <- comparison[,3]*1000000/sum(comparison[,3])
  
  mydata_count[,4] <- comparison[,4]*1000000/sum(as.numeric(comparison[,4]))
  mydata_count[,5] <- comparison[,5]*1000000/sum(as.numeric(comparison[,5]))
  
  temp <- mydata_count[mydata_count[,2] >= 1  & mydata_count[,3] >= 1,]
  
  head(temp)
  
  name_figure <- paste(i,".cpm.pdf", sep="")
  
  pdf(name_figure)

  plot(log10(temp[,2]),log10(temp[,3]),
       xlab=x[[1]][1],
       ylab=x[[1]][2],
       main="Counts Per Million (log10)", cex=.3, xlim=c(0,5), 
       ylim=c(0,5))
  my_reg <- lm(log10(temp[,3]) ~ log10(temp[,2]))
  points(log10(temp[,2]),log10(temp[,3]), col=ifelse(abs(log10(temp[,3]) - log10(temp[,2])) >= criteria_diff,'red', 'black'), cex=.3)
  abline(my_reg, col="blue", lwd=2)
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)
  
  
  text(3.5,2,paste("a=",coef_a,";","b=",coef_b), cex=0.8, col="blue")
  
  coef_sp <- cor(log10(temp[,3]), log10(temp[,2]), method="spearman")
  coef_ps <- cor(log10(temp[,3]), log10(temp[,2]), method="pearson")
  coef_kendall <- cor(log10(temp[,3]), log10(temp[,2]), method="kendall")
  
  text(3.5,1.8,paste("spearman =",round(coef_sp,4)), cex=0.8,col="blue")
  text(3.5,1.6,paste("pearson =",round(coef_sp,4)), cex=0.8, col="blue")
  text(3.5,1.4,paste("kendall =",round(coef_kendall,4)), cex=0.8, col="blue")
  text(3.5,1.2,paste("r² =",round(coef_sp*coef_sp,4)), cex=0.8, col="blue")
  
  dev.off()
  ################################################################
  ## Nucleotide per million - trimmed
  
  name_figure <- paste(i,".npm.trimmed.pdf", sep="")
  
  pdf(name_figure)
  
  temp <- mydata_count[mydata_count[,4] >= 1  & mydata_count[,5] >= 1,]
  
  plot(log10(temp[,4]),log10(temp[,5]),
       xlab=x[[1]][1],
       ylab=x[[1]][2],
       main="Nucleotides Per Million (log10) - trimmed reads", cex=.3, xlim=c(0,5), 
       ylim=c(0,5))
  my_reg <- lm(log10(temp[,5]) ~ log10(temp[,4]))
  points(log10(temp[,4]),log10(temp[,5]), 
         col=ifelse(abs(log10(temp[,5]) - log10(temp[,4])) >= criteria_diff,'red', 'black'), cex=.3)
  abline(my_reg, col="blue", lwd=2)
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)

  text(3.5,2,paste("a=",coef_a,";","b=",coef_b), cex=0.8, col="blue")
  
  coef_sp <- cor(log10(temp[,5]), log10(temp[,4]), method="spearman")
  coef_ps <- cor(log10(temp[,5]), log10(temp[,4]), method="pearson")
  coef_kendall <- cor(log10(temp[,5]), log10(temp[,4]), method="kendall")
  
  text(3.5,1.8,paste("spearman =",round(coef_sp,4)), cex=0.8,col="blue")
  text(3.5,1.6,paste("pearson =",round(coef_sp,4)), cex=0.8, col="blue")
  text(3.5,1.4,paste("kendall =",round(coef_kendall,4)), cex=0.8, col="blue")
  text(3.5,1.2,paste("r² =",round(coef_sp*coef_sp,4)), cex=0.8, col="blue")
  
  dev.off()
  
  ############################################################################
  #### Nucleotides per million - untrimmed
  
  ## Nucleotide per million
  
  ################################################################
  
  name_figure <- paste(i,".altman_blad.pdf", sep="")
  
  pdf(name_figure)
  
  plot(log10(rowMeans(temp[,2:3])),
       log10(temp[,3]) - log10(temp[,2]),
       xlab=paste("Mean CPM (log10)"),
    ylab=paste(x[[1]][2],"-",x[[1]][1],"(log10)"),
    ylim=c(-1.5,1.5),
    cex=0.3)
  
  my_mean <- mean(log10(temp[,3]) - log10(temp[,2]))
  my_sd <- sd(log10(temp[,3]) - log10(temp[,2]))
  
  abline(h=criteria_diff, col="red")
  abline(h=my_mean+my_sd, col="green")
  abline(h=my_mean, col="blue")
  abline(h=my_mean-my_sd, col="green")
  abline(h=(-1)*criteria_diff, col="red")

  dev.off()  

  # d2 <- mydata_count[log10(mydata_count$V7) - log10(mydata_count$V3) >= criteria_diff &
  #                      log10(mydata_count$V7) >= 1 & log10(mydata_count$V3) >= 1,]
  # d1 <- mydata_count[log10(mydata_count$V7) - log10(mydata_count$V3) <= (-1)*criteria_diff  &
  #                      log10(mydata_count$V7) >= 1 & log10(mydata_count$V3) >= 1,]
  # 
  # d4 <- mydata_count2[log10(mydata_count2$V8) - log10(mydata_count2$V4) >= criteria_diff &
  #                      log10(mydata_count2$V8) >= 1 & log10(mydata_count2$V4) >= 1,]
  # d3 <- mydata_count2[log10(mydata_count2$V8) - log10(mydata_count2$V4) <= (-1)*criteria_diff  &
  #                      log10(mydata_count2$V8) >= 1 & log10(mydata_count2$V4) >= 1,]
  # 
  # # print(paste("Mapped Reads = ",i, "-",x[[1]][1] ,":",nrow(d1),";",x[[1]][2] ,":", nrow(d2),";", "Total" ,":",nrow(d2) + nrow(d1), "---", raw_counts))
  # 
  # print(paste(i,"Reads vs Bases = ", nrow(d2)+nrow(d1),"vs", nrow(d3) + nrow(d4)))
  # 
  # 
  # if(nrow(d2) + nrow(d1) > 50){
  #   #boxplot(d1$V2,d2$V2,notch = T, outline = F, names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])))
  #   vioplot(log10(d1$V2), log10(d2$V2), names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])), col="white")
  #   beanplot(log10(d1$V2), log10(d2$V2), names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])))
  #   }
  # 
  # file_to_write <- paste(x[[1]][1],"_",x[[1]][2],".tx_comparison_3-fold", sep="")
  # 
  # if(file.exists(file_to_write)){
  #   command <- paste("rm",file_to_write)
  #   system(command)
  # }
  # 
  # if(nrow(d2) > 0){
  #   y1 <- cbind(d2[,1:2], x[[1]][2])
  #   write.table(y1,quote = F, file_to_write, append=T, sep = "\t", col.names = F, row.names = F)
  # }
  # 
  # if(nrow(d1) > 0){
  #   y2 <- cbind(d1[,1:2], x[[1]][1])
  #   write.table(y2,quote = F, file_to_write, append=T, sep = "\t", col.names = F , row.names = F)
  # }
}

##################################################################################
################################# NUCLEOTIDE COVERAGE

for(i in comparisons){
  x <- strsplit(i,"_");
  
  mydata_count <- comparison_list[[i]][comparison_list[[i]]$V4 > limit_c & comparison_list[[i]]$V8 > limit_c,]
  
  mydata_count$V4 <- mydata_count$V4*1000000/sum(as.numeric(mydata_count$V4))
  mydata_count$V8 <- mydata_count$V8*1000000/sum(as.numeric(mydata_count$V8))
  
  mydata_count <- mydata_count[mydata_count$V4 >= 1  & mydata_count$V8 >= 1,]
  
  raw_counts <- nrow(mydata_count)
  
  name_figure <- paste(i,".bpm.pdf", sep="")
  
  #pdf(name_figure)
  
  plot(log10(mydata_count$V4),log10(mydata_count$V8),
       xlab=x[[1]][1],
       ylab=x[[1]][2],
       main="Bases Per Million (log10)", cex=.3, xlim=c(0,5), 
       ylim=c(0,5))
  my_reg <- lm(log10(mydata_count$V8) ~ log10(mydata_count$V4))
  points(log10(mydata_count$V4),log10(mydata_count$V8), col=ifelse(abs(log10(mydata_count$V8) - log10(mydata_count$V4)) >= criteria_diff,'red', 'black'), cex=.3)
  abline(my_reg, col="blue", lwd=2)
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)
  
  
  text(3.5,2,paste("a=",coef_a,";","b=",coef_b), cex=0.8, col="blue")
  
  abline(my_reg, col="blue")
  
  coef_sp <- cor(log10(mydata_count$V8), log10(mydata_count$V4), method="spearman")
  coef_ps <- cor(log10(mydata_count$V8), log10(mydata_count$V4), method="pearson")
  coef_kendall <- cor(log10(mydata_count$V8), log10(mydata_count$V4), method="kendall")
  
  text(3.5,1.8,paste("spearman =",round(coef_sp,4)), cex=0.8,col="blue")
  text(3.5,1.6,paste("pearson =",round(coef_sp,4)), cex=0.8, col="blue")
  text(3.5,1.4,paste("kendall =",round(coef_kendall,4)), cex=0.8, col="blue")
  text(3.5,1.2,paste("r² =",round(coef_sp*coef_sp,4)), cex=0.8, col="blue")
  
  #dev.off()
  
  name_figure <- paste(i,".bpm.altman_blad.pdf", sep="")
  
  #pdf(name_figure)
  
  plot(log10(rowMeans(subset(mydata_count, select= c(V4,V8)))),
       log10(mydata_count$V8) - log10(mydata_count$V4),
       xlab=paste("Mean BPM (log10)"),
       ylab=paste(x[[1]][2],"-",x[[1]][1],"(log10)"),
       ylim=c(-1.5,1.5),
       cex=0.3)
  
  my_mean <- mean(log10(mydata_count$V8) - log10(mydata_count$V4))
  my_sd <- sd(log10(mydata_count$V8) - log10(mydata_count$V4))
  
  abline(h=criteria_diff, col="red")
  abline(h=my_mean+my_sd, col="green")
  abline(h=my_mean, col="blue")
  abline(h=my_mean-my_sd, col="green")
  abline(h=(-1)*criteria_diff, col="red")
  
  #dev.off()  
  
  d2 <- mydata_count[log10(mydata_count$V8) - log10(mydata_count$V4) >= criteria_diff &
                       log10(mydata_count$V8) >= 1 & log10(mydata_count$V4) >= 1,]
  d1 <- mydata_count[log10(mydata_count$V8) - log10(mydata_count$V4) <= (-1)*criteria_diff  &
                       log10(mydata_count$V8) >= 1 & log10(mydata_count$V4) >= 1,]
  
  print(paste("Mapped Bases = ",i, "-",x[[1]][1] ,":",nrow(d1),";",x[[1]][2] ,":", nrow(d2),";", "Total" ,":",nrow(d2) + nrow(d1), "---", raw_counts))
  
  if(nrow(d2) + nrow(d1) > 50){
    #boxplot(d1$V2,d2$V2,notch = T, outline = F, names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])))
    vioplot(log10(d1$V2), log10(d2$V2), names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])), col="white")
    beanplot(log10(d1$V2), log10(d2$V2), names=c(paste("Up in",x[[1]][1]),paste("Up in", x[[1]][2])))
  }
  
  file_to_write <- paste(x[[1]][1],"_",x[[1]][2],".tx_comparison_3-fold", sep="")
  
  if(file.exists(file_to_write)){
    command <- paste("rm",file_to_write)
    system(command)
  }
  
  if(nrow(d2) > 0){
    y1 <- cbind(d2[,1:2], x[[1]][2])
    write.table(y1,quote = F, file_to_write, append=T, sep = "\t", col.names = F, row.names = F)
  }
  
  if(nrow(d1) > 0){
    y2 <- cbind(d1[,1:2], x[[1]][1])
    write.table(y2,quote = F, file_to_write, append=T, sep = "\t", col.names = F , row.names = F)
  }
}


#####################################################################################
#### Normalize by size factors

library(bioconductor)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)

## read complete table

## calculate normalized factors

featureCounts_file <- read.table("allReplicates_output_allGenes", row.names = 1, header = T)

head(featureCounts_file)

cts <- featureCounts_file[,6:16]

colnames(cts) <- replicates

coldata <- read.table("colData", row.names = 1, header = T)

dds2 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds2

dds2 <- dds[ rowSums(counts(dds2)) > 1, ]

dds <- estimateSizeFactors(dds)
rld <- rlog(dds)

plotPCA(rld, intgroup=c("condition", "type"))

se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))

plotPCA( DESeqTransform( se ), intgroup=c("condition", "type"))

sampleDists <- dist(t(assay(rld)))

library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#################################################
# Plot nucleotide dist

for(i in comparisons){
  print(i)
  x <- strsplit(i,"_")
  
  limit_c <- 100;
  
  mydata_x <- comparison_list[[i]][comparison_list[[i]]$V3 > limit_c & comparison_list[[i]]$V7 > limit_c,]$V4
  mydata_y <- comparison_list[[i]][comparison_list[[i]]$V3 > limit_c & comparison_list[[i]]$V7 > limit_c,]$V8
  
  plot(log10(mydata_x),log10(mydata_y),
       xlab=x[[1]][1],
       ylab=x[[1]][2],
       main="Aligned Nucleotides", cex=.3)
  
  my_reg <- lm(log10(mydata_y) ~ log10(mydata_x))
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)
  
  
  text(1000000,500,paste("a=",coef_a,";","b=",coef_b), cex=0.8, col="red")
  
  abline(lm(log10(mydata_y) ~ log10(mydata_x)), col="red")
  
  coef_sp <- cor(log10(mydata_y), log10(mydata_x), method="spearman")
  coef_ps <- cor(log10(mydata_y), log10(mydata_x), method="pearson")
  coef_kendall <- cor(log10(mydata_y), log10(mydata_x), method="kendall")
  
  text(1000000,300,paste("spearman =",round(coef_sp,4)), cex=0.8)
  text(1000000,180,paste("pearson =",round(coef_sp,4)), cex=0.8)
  text(1000000,120,paste("kendall =",round(coef_kendall,4)), cex=0.8)
  
}


###################
###### calculate nt correlation

nt_correlation_file <- "nt_correlation"

if(file.exists(nt_correlation_file)){
  command <- paste("rm",nt_correlation_file)
  system(command)
}

for(i in comparisons){
  x <- strsplit(i,"_")
  
  limit_c <- 100;
  
  mydata_x <- comparison_list[[i]][comparison_list[[i]]$V3 > limit_c & comparison_list[[i]]$V7 > limit_c,]$V4
  mydata_y <- comparison_list[[i]][comparison_list[[i]]$V3 > limit_c & comparison_list[[i]]$V7 > limit_c,]$V8
  
  my_reg <- lm(log10(mydata_y) ~ log10(mydata_x))
  
  plot(log10(mydata_y) - log10(mydata_x), ylab=paste(x[[1]][2],"-",x[[1]][1]),cex=.3, col="lightgreen")
  abline(h=0,col="red")
  
  points(log10(mydata_y*1000000/sum(as.numeric(mydata_y))) - log10(mydata_x*1000000/sum(as.numeric(mydata_x))), ylab=paste(x[[1]][2],"-",x[[1]][1]),cex=.3, col="blue")
  
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)
  
  
  coef_sp <- cor(log10(mydata_y), log10(mydata_x), method="spearman")
  coef_ps <- cor(log10(mydata_y), log10(mydata_x), method="pearson")
  coef_kendall <- cor(log10(mydata_y), log10(mydata_x), method="kendall")
  
  line <- c(i, coef_a, coef_b, coef_sp, coef_ps, coef_kendall)
  
  write(line,nt_correlation_file, append = TRUE, sep = " ", ncolumns = length(line))
  
}



