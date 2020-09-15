## 
setwd("/home/vitor/PRJ.SRP064142/QUANTIFY_REPLICATES/")
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

file_terminal <- ".gene.length.count.total_nt.mean_nt.nt_by_length"

comparison_list <- list()

for(i in comparisons){
  name_comp <- paste(i,file_terminal,sep="")
  comparison_list[[i]] <- read.table(name_comp)
  
}

correlations <- list()




for(i in comparisons){
  x <- strsplit(i,"_")
  
  mydata_x <- comparison_list[[i]][comparison_list[[i]]$V3 > 0 & comparison_list[[i]]$V7 > 0,]$V3
  mydata_y <- comparison_list[[i]][comparison_list[[i]]$V3 > 0 & comparison_list[[i]]$V7 > 0,]$V7
  
  plot(mydata_x,mydata_y,
       xlab=x[[1]][1],
       ylab=x[[1]][2],
       log="xy",
       main="Read Counts", cex=.3)
  
  my_reg <- lm(log10(mydata_y) ~ log10(mydata_x))
  
  coef_a <- round(10^my_reg$coefficients[[2]][1],2)
  coef_b <- round(10^my_reg$coefficients[[1]][1],2)
  
  
  text(10600,80,paste("a=",coef_a,";","b=",coef_b), cex=0.8, col="red")
  
  abline(lm(log10(mydata_y) ~ log10(mydata_x)), col="red")
  
  coef_sp <- cor(log10(mydata_y), log10(mydata_x), method="spearman")
  coef_ps <- cor(log10(mydata_y), log10(mydata_x), method="pearson")
  coef_kendall <- cor(log10(mydata_y), log10(mydata_x), method="kendall")
  
  text(10600,50,paste("spearman =",round(coef_sp,4)), cex=0.8)
  text(10500,30,paste("pearson =",round(coef_ps,4)), cex=0.8)
  text(10500,20,paste("kendall =",round(coef_kendall,4)), cex=0.8)
  text(10500,10,paste("pearson² =",round(coef_ps^2,4)), cex=0.8)
  
  print(paste(i,coef_sp, coef_ps, coef_kendall, coef_ps^2))

}


