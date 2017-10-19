setwd("/home/vitor/Proj_ProC_R/mappings_star")

library("ggplot2")

count_tab <- read.table("count_table", header=T, row.names = 1)

head(count_tab)
num_replicates <- ncol(count_tab)
replicates_names <- colnames(count_tab)

for(i in 1:(num_replicates-1)){
  for(j in (i+1):num_replicates){
    
    rep1 <- replicates_names[i];
    rep2 <- replicates_names[j];
    
    my_comp <- cbind(count_tab[[rep1]],count_tab[[rep2]]);
    rownames(my_comp) <- rownames(count_tab);
    colnames(my_comp) <- c(rep1,rep2);
    
    cpms <- my_comp[my_comp[,1] >= 10 & my_comp[,2] >= 10,];
    
    cpms[,1] <- cpms[,1]*10^6/sum(cpms[,1]);
    
    cpms[,1] <- log2(cpms[,1]);
    
    cpms[,2] <- cpms[,2]*10^6/sum(cpms[,2]);
    
    cpms[,2] <- log2(cpms[,2]);
    
    
    output_file_name <- paste(rep1,"_",rep2,"cpms.png", sep="")
    png(output_file_name)
    
    cpms <- as.data.frame(cpms)
    
    p <- ggplot(cpms, aes(cpms[,1],cpms[,2])) 
    print(p + geom_point() + stat_smooth(method="lm", se=T, col="red") + labs(x=rep1,y=rep2))
    #+ annotate("text", label = "R² = 0.8987; p-value < 2.2e-16", x = 4, y = 1, size = 4, colour = "red")
    
    dev.off()
    
    }
}

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
