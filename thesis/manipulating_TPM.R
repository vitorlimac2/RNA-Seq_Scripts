workdir <- "/home/vitor"
setwd(workdir)

n <- read.table("TPM_normalized_frac.csv", header=FALSE, skip = 1, row.names = 1)

colnames(n) <- c("length","transcript","ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2","HiSeq1","Hiseq2")
r <- read.table("TPM_raw_frac.csv", header=FALSE, skip = 1, row.names = 1)
colnames(r) <- c("length","transcript","ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2","HiSeq1","Hiseq2")

top_200_r_proR1 <- rownames(head(r[order(r$ProR1,decreasing = T),], n=200))
top_200_n_proR1 <- rownames(head(n[order(n$ProR1,decreasing = T),], n=200))
top_200_n_Hiseq1 <-rownames(head(n[order(n$HiSeq1,decreasing = T),], n=200))



cor.test(top_200_n_proR1, top_200_r_proR1)


head(n[order(n$HiSeq1,decreasing = T),])
