
cor.my.ecdf <- function(REP1,REP2,STATUS = NULL, n_iterations, n_observations = NULL, norm_path, raw_path){
  
  
  vcor1 <- vector(mode="numeric", length=n_iterations)
  vcor2 <- vector(mode="numeric", length=n_iterations)
  
  status1 <- read.table(paste(REP1,".status",sep=""))
  
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
  
  keep <- status1[status1$V2 == STATUS,]$V1
  
  if(is.null(STATUS))
    y <- f1
  else
    y <- f1[keep,]
  i = 0;
  
  ## number of genes to sub-sampling
  
  if(is.null(n_observations))
    n_observations <- round(.3*nrow(y))
  while (i <= n_iterations) {
    f1_lf <-  y[sample(nrow(y), n_observations), ]
    x <- merge(f1_lf, f2, by = "gene", sort = T)
    #  print(nrow(x))
    c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
    c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
    
    vcor1[i] <- c1$estimate^2
    vcor2[i] <- c2$estimate^2
    i <- i + 1
  }
  
  
  return(t.test(vcor2,vcor1,paired=TRUE,alternative="greater")$p.value)
  
}

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")


r <- 50
n <- r*6
vcorC1 <- vector(mode="numeric", length=n)
vcorC2 <- vector(mode="numeric", length=n)
vcorC3 <- vector(mode="numeric", length=n)
vcorC4 <- vector(mode="numeric", length=n)
vcorC5 <- vector(mode="numeric", length=n)
vcorR1 <- vector(mode="numeric", length=n)
vcorR2 <- vector(mode="numeric", length=n)

replicates <- c("ProC1", "ProC2", "ProC3", "ProC4",
                "ProC5", "ProR1", "ProR2", "HiSeq1", "HiSeq2")


########## ProC vs ProR
## for cada replica
for(j in 1:5){
  z <- 1
  for(i in 1:r){
    
    for(k in 8:9){
      REP1 <- replicates[j]
      REP2 <- replicates[k]
      if(replicates[j]==replicates[k])
        next;
      if(replicates[j]=="HiSeq1" | replicates[j]=="HiSeq2")
        next;
      if(replicates[j]=="ProC1"){
        vcorC1[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
      }
      if(replicates[j]=="ProC2"){
        vcorC2[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
        
      }
      if(replicates[j]=="ProC3"){
        vcorC3[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 200, norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
        
      }
      if(replicates[j]=="ProC4"){
        vcorC4[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
      }
      if(replicates[j]=="ProC5"){
        vcorC5[5] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
        
      }
      if(replicates[j]=="ProR1"){
        vcorR1[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
        
      }
      if(replicates[j]=="ProR2"){
        vcorR2[z] <- cor.my.ecdf(REP1, REP2, 
                                 n_iterations = 100, 
                                 n_observations = 200,
                                 norm_path = "NormalizedAll", 
                                 raw_path = "RawAll")
        
      }
      z <- z + 1
    }
  }
}


#### APENAS PLOT

qvcorC1 <- p.adjust(vcorC1, method = "fdr",n=length(vcorC1))
qvcorC2 <- p.adjust(vcorC2, method = "fdr",n=length(vcorC2))
qvcorC3 <- p.adjust(vcorC3, method = "fdr",n=length(vcorC3))
qvcorC4 <- p.adjust(vcorC4, method = "fdr",n=length(vcorC4))
qvcorC5<- p.adjust(vcorC5, method = "fdr",n=length(vcorC5))

plot_range <- range(c(qvcorC1,qvcorC2,qvcorC3,qvcorC4,qvcorC5))

plot(main = "Fragmentação Química vs Enzimática",
     ecdf(qvcorC1), 
     lwd=2, 
     col="black",
     xlab="FDR",
     ylab="Sub-amostragem (%)",
     yaxt='n')
axis(2, at= c(0,0.2,.4,.60,.80,1),labels=c(0,20,40,60,80,100), col.axis="black", las=2)

plot(main = "Fragmentação Química vs Enzimática",
     ecdf(qvcorC6), 
     lwd=2, 
     col="red",
     xlab="FDR",
     ylab="Sub-amostragem (%)",
     yaxt='n')
axis(2, at= c(0,0.2,.4,.60,.80,1),labels=c(0,20,40,60,80,100), col.axis="black", las=2)

qvcorR1 <- p.adjust(vcorR1, method = "fdr",n=length(vcorR1))
qvcorR2<- p.adjust(vcorR2, method = "fdr",n=length(vcorR2))

max(qvcorR1, qvcorR2)

plot(ecdf(qvcorC2),  lwd=2, add=TRUE, lty="dashed", col="black")
plot(ecdf(qvcorC3),  lwd=2, add=TRUE, lty="dashed", col="blue")
plot(ecdf(qvcorC4),  lwd=2, add=TRUE, lty="dashed", col="darkgreen")
plot(ecdf(qvcorC5),  lwd=2, add=TRUE, lty="dashed", col="darkyellow")
