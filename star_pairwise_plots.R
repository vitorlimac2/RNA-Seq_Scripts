setwd("/home/vitor/Proj_ProC_R/mappings_star")

library("ggplot2")
library("DESeq2")

##########################################################################################
#### Gene Expression with Reads

count_tab <- read.table("count_table", header=T, row.names = 1)

head(count_tab)
num_replicates <- ncol(count_tab)
replicates_names <- colnames(count_tab)


read_mean_hash <- new.env(hash=T, parent=emptyenv())
mean_read_length <- c(153,138,155,136,155,137,158,140,160,94,92)

for(i in 1:length(replicates_names)){
  read_mean_hash[[replicates_names[i]]]<-mean_read_length[i]
  print(read_mean_hash[[replicates_names[i]]])
}




for(i in 1:(num_replicates-1)){
  for(j in (i+1):num_replicates){
    
    rep1 <- replicates_names[i];
    rep2 <- replicates_names[j];
    
    my_comp <- cbind(count_tab[[rep1]],count_tab[[rep2]]);
    rownames(my_comp) <- rownames(count_tab);
    colnames(my_comp) <- c(rep1,rep2);
    
    cpms <- my_comp[my_comp[,1] >= 1 & my_comp[,2] >= 1,];
    
    cpms[,1] <- cpms[,1]*10^9*(read_mean_hash[[rep1]]/read_mean_hash[[rep2]])/sum(cpms[,1]);
    
    cpms[,1] <- log2(cpms[,1]);
    
    cpms[,2] <- cpms[,2]*10^9*(read_mean_hash[[rep2]]/read_mean_hash[[rep1]])/sum(cpms[,2]);
    
    cpms[,2] <- log2(cpms[,2]);
    
    
    output_file_name <- paste(rep1,"_",rep2,"normalized_by_r1_under_r2.png", sep="")
    png(output_file_name)
    
    cpms <- as.data.frame(cpms)
    
    p <- ggplot(cpms, aes(cpms[,1],cpms[,2])) 
    print(p + geom_point() + stat_smooth(method="lm", se=T, col="red") + labs(x=rep1,y=rep2))
    #+ annotate("text", label = "R² = 0.8987; p-value < 2.2e-16", x = 4, y = 1, size = 4, colour = "red")
    
    dev.off()
    
    }
}


##########################################################################################
#### Gene Expression with Nucleotides

nt_tab <- read.table("nucleotides_table.csv", header=T, row.names = 1)

head(nt_tab)
num_replicates <- ncol(nt_tab)
replicates_names <- colnames(nt_tab)

for(i in 1:(num_replicates-1)){
  for(j in (i+1):num_replicates){
    
    rep1 <- replicates_names[i];
    rep2 <- replicates_names[j];
    
    my_comp <- cbind(nt_tab[[rep1]],nt_tab[[rep2]]);
    rownames(my_comp) <- rownames(nt_tab);
    colnames(my_comp) <- c(rep1,rep2);
    
    head(my_comp)
    
    nms <- my_comp[my_comp[,1] >= 10 & my_comp[,2] >= 10,];
    
    nms[,1] <- nms[,1]*10^6/sum(as.numeric(nms[,1]));
    
    nms[,1] <- log2(nms[,1]);
    
    nms[,2] <- nms[,2]*10^6/sum(as.numeric(nms[,2]));
    
    nms[,2] <- log2(nms[,2]);
    
    
    output_file_name <- paste(rep1,"_",rep2,"nms.png", sep="")
    png(output_file_name)
    
    nms <- as.data.frame(nms)
    
    p <- ggplot(nms, aes(nms[,1],nms[,2])) 
    print(p + geom_point() + stat_smooth(method="lm", se=T, col="red") + labs(x=rep1,y=rep2))
    #+ annotate("text", label = "R² = 0.8987; p-value < 2.2e-16", x = 4, y = 1, size = 4, colour = "red")
    
    dev.off()
    
  }
}

######################################################################################
##################### DESeq2

cts <- count_tab

keep = rowSums(count_tab > 1) >= 10

cts = count_tab[keep,]

head(cts)

condition <- c("C","C","C","C","C","C","C","C","C", "R","R")
initialRNA <- c(2,2,2,2,2,0.2,0.2,0.2,0.2,2,2)
colData <- data.frame(row.names = replicates_names, condition,initialRNA)
colData

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ initialRNA + condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=c("C","R"))
dds <- DESeq(dds)

rld <- rlogTransformation(dds, blind=FALSE)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE)),
                           colData=colData(dds))
plotPCA(DESeqTransform(se), intgroup=c("condition","initialRNA")) + 
  geom_text(aes(label = replicates_names), position = position_nudge(y = 1))
plotPCA(rld, intgroup=c("condition","initialRNA")) + 
  geom_text(aes(label = replicates_names), position = position_nudge(y = 1))


###########################################################################################
#### Mapping rate per Read Length

mapping_summary_dir <- "/home/vitor/Proj_ProC_R/mappings_star/"

mapping_summary_file_terminal <- "Aligned.out.bam.read_map.status_length.sorted.filtered"

mapping_summary_files <- paste(replicates_names,mapping_summary_file_terminal,sep="")


i <- 0
for(my_file_name in mapping_summary_files){
  
  i <- i + 1
  
  file_path <- paste(mapping_summary_dir,my_file_name,sep="")
  my_tb <- read.table(file_path)
  print(file_path)
  # 1st plot with read length
  l_50 <- nrow(my_tb[my_tb[,2] < 50,])
  l_75 <- nrow(my_tb[my_tb[,2] >= 50 & my_tb[,2] < 75,])
  l_100 <- nrow(my_tb[my_tb[,2] >= 75 & my_tb[,2] < 100,])
  l_125 <- nrow(my_tb[my_tb[,2] >= 100 & my_tb[,2] < 125,])
  l_150 <- nrow(my_tb[my_tb[,2] >= 125 & my_tb[,2] < 150,])
  l_175 <- nrow(my_tb[my_tb[,2] >= 150 & my_tb[,2] < 175,])
  l_200 <- nrow(my_tb[my_tb[,2] >= 175 & my_tb[,2] < 200,])
  l_225 <- nrow(my_tb[my_tb[,2] >= 225 & my_tb[,2] < 250,])
  l_250 <- nrow(my_tb[my_tb[,2] >= 250,])
  
  l_total <- nrow(my_tb)
  
  m_50 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] < 50,])
  m_75 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 50 & my_tb[,2] < 75,])
  m_100 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 75 & my_tb[,2] < 100,])
  m_125 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 100 & my_tb[,2] < 125,])
  m_150 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 125 & my_tb[,2] < 150,])
  m_175 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 150 & my_tb[,2] < 175,])
  m_200 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 175 & my_tb[,2] < 200,])
  m_225 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 225 & my_tb[,2] < 250,])
  m_250 <- nrow(my_tb[my_tb[,1]==1 & my_tb[,2] >= 250,])
  
  l_vector <- c(l_50, l_75, l_100, l_125, l_150, l_175, l_200, l_225, l_250)*100/l_total
  
  m_total <- nrow(my_tb[my_tb[,1]>=1, ])
  
  m_vector <- c(m_50, m_75, m_100, m_125, m_150, m_175, m_200, m_225, m_250, m_275)*100/m_total
  
  intervals <- c("50","75","100","125","150","175","200","225",">=250")

  png(paste(replicates_names[i],"total_vs_mapped.png", sep=""))
  plot(main=replicates_names[i], m_vector, type="l",col="blue", axes=F, ylab="%", 
       xlab="Read length (bp)")
  
  legend(title="Reads","topright",c("Total","Unique Mapped"), fill=c("black","blue"))

  axis(1,lwd=2,labels=intervals, at=1:9)
  axis(2,lwd=2)
  lines(l_vector)
  
  dev.off()
  
  #png(paste(main_title,"diff_map.png", sep=""))
  
  #plot(main=main_title,m_vector-l_vector, axes=F, type="l",xlab="Read Length", ylab="(% Mapped) - (% Total)")
  #axis(1,lwd=2,labels=intervals, at=1:10)
  #abline(a=0,b=0, col="darkred")
  #axis(2,lwd=2)
  
  #dev.off()
  
  # 2nd plot with mapping rate per length
  
  rm(my_tb)
}

