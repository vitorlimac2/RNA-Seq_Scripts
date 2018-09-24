workdir <- "/home/vitor/Downloads"

setwd(workdir)

gene1 <- read.table("ENSG00000255823.HiSeq1.InsertSizeDistribution")$V1
gene2 <- read.table("ENSG00000269028.HiSeq1.InsertSizeDistribution")$V1

rep <- read.table("HiSeq1.InsertSizeDistribution")$V1
my_breaks <- pretty(range(c(gene1,gene2,rep)), n = 20)

h1 <- hist(gene1, breaks=my_breaks, plot=FALSE)
h2 <- hist(gene2, breaks=my_breaks, plot=FALSE)
hrep <- hist(rep, breaks=my_breaks, plot=FALSE)

max_y <- max(h1$counts*100/length(gene1), 
             h2$counts*100/length(gene2),
             hrep$counts*100/length(rep))

plot(main = "Sample HiSeq1", 
     hrep$mids, hrep$counts*100/length(rep), 
     type="l", 
     ylab = "%", 
     xlab = "Insert Size (bp)", 
     col="black", 
     ylim=c(0,max_y), xlim=c(0,350), lwd=2)

lines(h2$mids, h1$counts*100/length(gene1), col="red", lty=1, lwd=2)
lines(h1$mids, h1$counts*100/length(gene1), col="blue", lty=2, lwd=2)

legend("topright", 
       legend = c("HiSeq1", "ENSG00000255823", "ENSG00000269028"), 
       col = c("black", "blue","red"), bty = 'n', cex=1, lty = c(1,2,2))
