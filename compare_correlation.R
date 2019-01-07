## read file

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")

f1 <- read.table("ProR1.Analisys")
f2 <- read.table("ProC3.Analisys")

y <- f1[f1$V4=="L",]

n_iterations <- 500;
i = 0;

v1 <- vector(mode="numeric", length=n_iterations)
v2 <- vector(mode="numeric", length=n_iterations)

while (i <= n_iterations) {
  f1_lf <-  y[sample(nrow(y), 200), ]
  x <- merge(f1_lf, f2, by = c("V1","V1"))
  print(nrow(x))
  c1 <- cor.test(x$V2.x, x$V2.y, method="spearman")
  c2 <- cor.test(x$V3.x, x$V3.y, method="spearman")
  if(c1$p.value <= 0.01 &  c1$p.value <= 0.01){
    v1[i] <- c1$estimate^2
    v2[i] <- c2$estimate^2
    i <- i + 1
  }else{
    i <- i -1
  }
}


plot(main = "Distribuição Empírica Cumulativa - ProR2 vs ProC1 - Baixa Fragmentação",
     ecdf(v1), lwd=2,xlim=range(c(v1, v2)), col="red",
     xlab="Correlação (Spearman)")
plot(ecdf(v2),  lwd=2, add=TRUE, lty="dashed", col="blue")

legend(min(v1,v2), .9,box.lwd = 0,box.col = "white",bg = "white", 
       legend = c("Não-normalizado", "Normalizado" ), 
       fill=c("red", "blue"))

####################################################################3
