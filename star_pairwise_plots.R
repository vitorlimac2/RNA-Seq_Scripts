setwd("/home/vitor/Proj_ProC_R/mappings_star")

count_tab <- read.table("count_table", header=T, row.names = 1)

head(count_tab)
num_replicates <- ncol(count_tab)
replicates_names <- colnames(count_tab)

for(i in 1:(num_replicates-1)){
  for(j in (i+1):num_replicates){
    print(paste(i,j, sep = " "));
    
    rep1 <- replicates_names[i];
    rep2 <- replicates_names[j];
    
    my_comp <- cbind(count_tab[[rep1]],count_tab[[rep2]])
    rownames(my_comp) <- rownames(count_tab)
    colnames(my_comp) <- c(rep1,rep2)
    
    head(my_comp)
    
    cpms <- my_comp[my_comp[,1] != 0 & my_comp[,2] != 0,]
    
    head(cpms)
    
    }
}
