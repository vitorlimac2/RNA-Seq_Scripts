# ## read file
# 
# ### Primeiro precisa testar se a distribuição das amostras é bem diferente
# ### com KS
# 
# mse <- function(vec1,vec2,x){
#   lm1 <- lm(vec1~vec2,x,);
#   sm <- summary(lm1);
#   return(mean(sm$residuals^2))
#   }
# 
# setwd("/media/vitor/Seagate Expansion Drive/Thesis/")
# 
# 
# f1 <- read.table("ProC1.Analisys")
# f2 <- read.table("ProC2.Analisys")
# 
# f1$V2 <- f1$V2*1000000/sum(f1$V2)
# f1$V3 <- f1$V3*1000000/sum(f1$V3)
# f1
# f2$V2 <- f2$V2*1000000/sum(f2$V2)
# f2$V3 <- f2$V3*1000000/sum(f2$V3)
# 
# f1$V2 <- log(f1$V2 + 0.1)
# f1$V3 <- log(f1$V3 + 0.1)
# f2$V2 <- log(f2$V2 + 0.1)
# f2$V3 <- log(f2$V3 + 0.1)
# 
# y <- f1[f1$V4 == "N",]
# y
# 
# n_iterations <- 100;
# i = 0;
# 
# v1 <- vector(mode="numeric", length=n_iterations)
# v2 <- vector(mode="numeric", length=n_iterations)
# mse1 <- vector(mode="numeric", length=n_iterations)
# mse2 <- vector(mode="numeric", length=n_iterations)
# 
# while (i <= n_iterations) {
#   f1_lf <-  y[sample(nrow(y), 200), ]
#   x <- merge(f1_lf, f2, by = c("V1","V1"))
#   print(nrow(x))
#   c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
#   c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
#   if(c1$p.value <= 0.01 &  c1$p.value <= 0.01){
#     mse1[i] <- mse(x$V2.x,x$V2.y,x)
#     print(mse(x$V2.x,x$V2.y,x))
#     mse2[i] <- mse(x$V3.x, x$V3.y,x)
#     v1[i] <- c1$estimate^2
#     v2[i] <- c2$estimate^2
#     i <- i + 1
#   }else{
#     i <- i -1
#   }
# }
# 
# 
# plot(main = "Distribuição Empírica Cumulativa - ProR2 vs ProC1 - Baixa Fragmentação",
#      ecdf(v1), lwd=2,xlim=range(c(v1, v2)), col="red",
#      xlab="Correlação (Spearman)",
#      ylab="%")
# plot(ecdf(v2),  lwd=2, add=TRUE, lty="dashed", col="black")
# 
# legend(title = "CPM (log)",min(v1,v2), .9,box.lwd = 0,box.col = "white",bg = "white", 
#        legend = c("Não-ajustado", "Ajustado" ), 
#        fill=c("red", "black"))
# 
# plot(ylab="Erro Quadrático Médio", 
#      col='red',
#      log='y',mse1,ylim=range(c(mse1, mse2)), 
#      xlab='Reamostragem',
#      type='l', lwd=3)
# lines(mse2, col='black',type='l', lwd=3)
# boxplot(mse1,mse2, outline = F)



####################################################################3

mse <- function(vec1,vec2,x){
   lm1 <- lm(vec1~vec2,x,);
   sm <- summary(lm1);
   return(mean(sm$residuals^2))
   }


plot.my.ecdf <- function(REP1,REP2,STATUS, n_iterations, n_observations, norm_path, raw_path){
 
  ## define replicate to compare
  #REP1 <- "ProR1"
  #REP2 <- "ProC3"
  ## read status file
  
  status1 <- read.table(paste(REP1,".status",sep=""))
  
  #status2 <- read.table(paste(REP2,".status",sep=""))

  ## read raw file
  raw_file <- read.table(raw_path, 
                         header = T)
  ## read norm file
  norm_file <- read.table(norm_path,
                          header=T)
  
  headings <- c("gene","ProC1", "ProC2", "ProC3", "ProC4",
                "ProC5", "ProC6", "ProC7", "ProC8",
                "ProC9", "ProR1", "ProR2", "HiSeq1", "HiSeq2")
  
  colnames(raw_file) <- headings
  colnames(norm_file) <- headings
  
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
  
  keep <- status1[status1$V2 == STATUS,]$V1
  y <- f1[keep,]
  y <- f1
  i = 0;
  
  vcor1 <- vector(mode="numeric", length=n_iterations)
  vcor2 <- vector(mode="numeric", length=n_iterations)
  mse1 <- vector(mode="numeric", length=n_iterations)
  mse2 <- vector(mode="numeric", length=n_iterations)
  variation1 <- vector(mode="numeric", length=n_iterations)
  variation2 <- vector(mode="numeric", length=n_iterations)
  
  while (i <= n_iterations) {
    f1_lf <-  y[sample(nrow(y), n_observations), ]
    x <- merge(f1_lf, f2, by = "gene", sort = T)
    print(nrow(x))
    c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
    c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
    if(c1$p.value <= 0.01 &  c1$p.value <= 0.01){
      mse1[i] <- mse(x$V2.x,x$V2.y,x)
      print(mse(x$V2.x,x$V2.y,x))
      mse2[i] <- mse(x$V3.x, x$V3.y,x)
      vcor1[i] <- c1$estimate^2
      vcor2[i] <- c2$estimate^2
      
      i <- i + 1
    }else{
      i <- i -1
    }
  }
  
  fragmentation_status = "Sem perfil definido";
  if(STATUS=="L"){
    fragmentation_status = "Baixa Fragmentação"
    
  }
  if(STATUS=="H"){
    fragmentation_status = "Alta Fragmentação"
  }
  
  plottitle <- paste(REP1," (",fragmentation_status,")"," vs ", REP2, sep="")
  
  plot(main = plottitle,
       ecdf(vcor1), 
       lwd=2, 
       xlim=range(c(vcor1, vcor2)), 
       col="red",
       xlab="Correlação (Spearman)",
       ylab="Sub-amostragem (%)",
       yaxt='n')
  
  axis(2, at= c(0,0.2,.4,.60,.80,1),labels=c(0,20,40,60,80,100), col.axis="black", las=2)
  
  print(ks.test(vcor1,vcor2, alternative = "two.sided"))
  
  plot(ecdf(vcor2),  lwd=2, add=TRUE, lty="dashed", col="black")
  
  legend(title = "CPM (log)",min(vcor1,vcor2), .9,box.lwd = 0,box.col = "white",bg = "white", 
         legend = c("Não-ajustado", "Ajustado" ), 
         fill=c("red", "black"), cex = 1)
  
   plot(ylab="Erro Quadrático Médio (CPM - log)", 
        col='red',
        log='y',mse1,ylim=range(c(mse1, mse2)), 
        xlab='Sub-amostragem',
        type='l', lwd=3)
   lines(mse2, col='black',type='l', lwd=3)
   
   boxplot(variation1,variation2, outline = F)

}

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")

plot.my.ecdf(REP1 = "ProC1", REP2 = "HiSeq1", 
             n_iterations = 100, 
             STATUS = "N", 
             n_observations = 200, 
             norm_path = "NormalizedAll", 
             raw_path = "RawAll")

library(debrowser)
startDEBrowser()
