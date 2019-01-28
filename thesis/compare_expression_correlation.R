
plot.my.ecdf <- function(REP1,REP2,STATUS = NULL, n_iterations, n_observations = NULL, norm_path, raw_path){
  
  
  vcor1 <- vector(mode="numeric", length=n_iterations)
  vcor2 <- vector(mode="numeric", length=n_iterations)
  
  #status1 <- read.table(paste(REP1,".status",sep=""))
  
  ## read raw file
  raw_file <- read.table(raw_path, 
                         header = T)
  ## read norm file
  norm_file <- read.table(norm_path,
                          header=T)
  
  m <- merge(raw_file, norm_file, by= c("gene","gene"), sort = T)
  
  f1 <- data.frame(gene = m$gene, V2 = m[,paste(REP1,".x",sep="")], 
                   V3= m[,paste(REP1,".y",sep="")])
  f2 <- data.frame(gene = m$gene, V2 = m[,paste(REP2,".x",sep="")], 
                   V3= m[,paste(REP2,".y",sep="")])
  
  f1$V2 <- f1$V2*1000000/sum(f1$V2)
  f1$V3 <- f1$V3*1000000/sum(f1$V3)
  #f1
  f2$V2 <- f2$V2*1000000/sum(f2$V2)
  f2$V3 <- f2$V3*1000000/sum(f2$V3)
  #f2
  f1$V2 <- log(f1$V2 + 0.1)
  f1$V3 <- log(f1$V3 + 0.1)
  f2$V2 <- log(f2$V2 + 0.1)
  f2$V3 <- log(f2$V3 + 0.1)
  
  
  #STATUS = "N"
  
  y <- f1
  i = 0;
  
  ## number of genes to sub-sampling
  
  if(is.null(n_observations))
    n_observations <- round(.25*nrow(y))
  
  print(paste("total =", nrow(y),";","sub-amostras=",n_observations))
  
  while (i <= n_iterations) {
    f1_lf <-  y[sample(nrow(y), n_observations), ]
    x <- merge(f1_lf, f2, by = "gene", sort = T)
    c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
    c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
    
    vcor1[i] <- c1$estimate^2
    vcor2[i] <- c2$estimate^2
    i <- i + 1
    # }else{
    #    i <- i - 1
    #  }
  }
  
  print(t.test(vcor2,vcor1,paired = TRUE, alternative="greater"))
  
  fragmentation_status = "Sem perfil definido";
  if(is.null(STATUS))
    fragmentation_status = ""
  else if(STATUS=="L"){
    fragmentation_status = "(Baixa Fragmentação)"
    
  }else if(STATUS=="H"){
    fragmentation_status = "(Alta Fragmentação)"
  }
  
  plottitle <- paste(REP1,fragmentation_status," vs ", REP2, sep="")
  
  plot(main = plottitle,
       ecdf(vcor1), 
       lwd=2, 
       xlim=range(c(vcor1, vcor2)), 
       col="red",
       xlab="Correlação (Spearman)",
       ylab="Sub-amostragem (%)",
       yaxt='n')
  
  axis(2, at= c(0,0.2,.4,.60,.80,1),labels=c(0,20,40,60,80,100), col.axis="black", las=2)
  
  plot(ecdf(vcor2),  lwd=2, add=TRUE, lty="dashed", col="black")
  
  legend(title = "CPM (log)",min(vcor1,vcor2), .9,box.lwd = 0,box.col = "white",bg = "transparent", 
         legend = c("Não-ajustado", "Ajustado" ), 
         fill=c("red", "black"), cex = 1)
  
  # plot(ylab="Erro Quadrático Médio (CPM - log)", 
  #      col='red',
  #      log='y',mse1,ylim=range(c(mse1, mse2)), 
  #      xlab='Sub-amostragem',
  #      type='l', lwd=3)
  # lines(mse2, col='black',type='l', lwd=3)
  
  return(t.test(vcor2,vcor1,paired=TRUE,alternative="greater"))
  
}

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")

#### APENAS PLOT
plot.my.ecdf(REP1 = "ProC", REP2 = "HiSeq", 
             n_iterations = 300, 
             #      STATUS = "L", 
             #       n_observations = 100, 
             norm_path = "NormalizedAllMean.tsv", 
             raw_path = "RawAllMean.tsv")
