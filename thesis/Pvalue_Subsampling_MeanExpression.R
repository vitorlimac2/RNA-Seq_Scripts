
cor.my.ecdf <- function(REP1,REP2,
                        n_iterations, 
                        n_observations = NULL, 
                        norm_path, 
                        raw_path){
  
  vcor1 <- vector(mode="numeric", length=n_iterations)
  vcor2 <- vector(mode="numeric", length=n_iterations)
  
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
  
  ## number of genes to sub-sampling
  
  y <- f1
  
  if(is.null(n_observations))
    n_observations <- round(.25*nrow(y))
  while (i <= n_iterations) {
    f1_lf <-  y[sample(nrow(y), n_observations), ]
    x <- merge(f1_lf, f2, by = "gene", sort = T)
    #  print(nrow(x))
    c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
    c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
    vcor1[i] <- c1$estimate
    vcor2[i] <- c2$estimate
    i <- i + 1
  }
  return(wilcox.test(vcor2,vcor1,paired = TRUE, alternative="greater")$p.value)
}

multiple.test <- function(REP1, REP2, norm_path, raw_path, n_iterations_cor, n_iterations){
  i <- 1;
  vec_Pvalues <- vector(mode="numeric", length = n_iterations)
  
  while(i <= n_iterations){
    vec_Pvalues[i] <- cor.my.ecdf(REP1 = REP1, 
                                  REP2 = REP2, 
                                n_iterations = n_iterations_cor,
                                norm_path = norm_path, 
                                raw_path = raw_path)
    i <- i + 1
  }
  return(vec_Qvalues <- p.adjust(vec_Pvalues, method = "fdr",n=length(vec_Pvalues)))
} 

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")


#replicates <- c("ProC_2",	"ProR",	"HiSeq")

qvcorC_2_HiSeq <- multiple.test(REP1 = "ProC_2", 
                                REP2 = "HiSeq",
                                n_iterations_cor = 100,
                                n_iterations = 200,
                                norm_path = "NormalizedAllMean.tsv",
                                raw_path = "RawAllMean.tsv")

qvcorC_2_R <- multiple.test(REP1 = "ProC_2", 
                            REP2 = "ProR",
                            n_iterations_cor = 100,
                            n_iterarions = 200,
                            norm_path = "NormalizedAllMean.tsv",
                            raw_path = "RawAllMean.tsv")

qvcorR_HiSeq <- multiple.test(REP1 = "ProR", 
                              REP2 = "HiSeq",
                              n_iterations_cor = 100,
                              n_iterarions = 200,
                              norm_path = "NormalizedAllMean.tsv",
                              raw_path = "RawAllMean.tsv")

my_breaks <- pretty(range(c(-log10(qvcorR_HiSeq), -log10(qvcorC_2_R), -log10(qvcorC_2_HiSeq))), n=30)

h_1 <- hist(-log10(qvcorC_2_HiSeq), breaks=my_breaks, plot=FALSE)
h_2 <- hist(-log10(qvcorC_2_R), breaks=my_breaks, plot=FALSE)
h_3 <- hist(-log10(qvcorR_HiSeq), breaks=my_breaks, plot=FALSE)


max_y <- max(h_1$counts*100/length(qvcorC_2_HiSeq), 
             h_2$counts*100/length(qvcorC_2_R),
             h_3$counts*100/length(qvcorR_HiSeq))

plot(main = "", 
     h_1$mids, 
     h_1$counts*100/length(qvcorC_2_HiSeq), 
     type="l", 
     ylab = "%", 
     xlab = "P-valor ajustado (FDR)", 
     col="black", 
     ylim=c(0,max_y), 
     cex.main=0.8,
     log="x",
     xaxt="n"
     )

lines(h_2$mids, h_2$counts*100/length(qvcorC_2_R), col="red")
lines(h_3$mids, h_3$counts*100/length(qvcorR_HiSeq), col="blue")

xaxis <- c(0.05,0.1,.2,.5,1,2)


axis(1, at=xaxis, labels = 10^(-1)*xaxis, las=2)


## Plot P-value adjusted #######################################
my_breaks <- pretty(-log10(qvcorC_2_R), 20)
hist(main = "ProC_2 vs HiSeq \n Aumento da correlação com ajuste",
     -log10(qvcorC_2), 
     breaks = my_breaks,
     xlab="P-valor ajustado (FDR)",
     ylab = "Frequência",
     xaxt="n")
axis(1, at=my_breaks, labels = 10^(-1)*my_breaks, las=2)

################################################################
