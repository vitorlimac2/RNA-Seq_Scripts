workdir <- "/home/vitor"
setwd(workdir)

n <- read.table("TPM_normalized_frac.csv", header=FALSE, skip = 1, row.names = 1)

colnames(n) <- c("length","transcript","ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2","HiSeq1","HiSeq2")
r <- read.table("TPM_raw_frac.csv", header=FALSE, skip = 1, row.names = 1)
colnames(r) <- c("length","transcript","ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2","HiSeq1","HiSeq2")

proReplicates <- c("ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2")

for(s in proReplicates){
  replicate <- s
  y = 50;
  
  n[,replicate]
  
  x <-  data.frame(n = n[n[,replicate] >= y & r[,replicate] >= y & n$HiSeq1 >= y & n$HiSeq2 >= y, replicate],
                   r = r[n[,replicate] >= y & r[,replicate] >= y & n$HiSeq1 >= y & n$HiSeq2 >= y,replicate],
                   h1= n[n[,replicate] >= y & r[,replicate] >= y & n$HiSeq1 >= y & n$HiSeq2 >= y,]$HiSeq1,
                   h2= n[n[,replicate] >= y & r[,replicate] >= y & n$HiSeq1 >= y & n$HiSeq2 >= y,]$HiSeq2)
  rownames(x) <- rownames(n[n[,replicate] >= y & r[,replicate] >= y & n$HiSeq1 >= y & n$HiSeq2 >= y,])
  
  x$h1_error_n <- ifelse(x$h1 > x$n, (x$h1 - x$n)*100/x$h1,(x$n - x$h1)*100/x$n)
  x$h1_error_r <- ifelse(x$h1 > x$r, (x$h1 - x$r)*100/x$h1,(x$r - x$h1)*100/x$r)
  x$h2_error_n <- ifelse(x$h2 > x$n, (x$h2 - x$n)*100/x$h2,(x$n - x$h2)*100/x$n) 
  x$h2_error_r <- ifelse(x$h2 > x$r, (x$h2 - x$r)*100/x$h2,(x$r - x$h2)*100/x$r)
  
  
  my_breaks <- pretty(range(c(x$h1_error_n, x$h1_error_r)), n=50)
  
  hist_n <- hist(x$h1_error_n, breaks=my_breaks, plot=FALSE)
  hist_r <- hist(x$h1_error_r, breaks=my_breaks, plot=FALSE)
  
  #max_y <- max(hist_n$counts*100/length(x$h1_error_n), hist_r$counts*100/length(x$h1_error_r))
  max_y <- 6
  
  png(paste(replicate,"png",sep="."))
  plot(main = replicate, hist_n$mids, hist_n$counts*100/length(x$h1_error_n), type="l", ylab = "%", xlab = "Error (%)", col="blue", ylim=c(0,max_y)) 
  
  lines(hist_r$mids, hist_r$counts*100/length(x$h1_error_r), col="red")
  
  legend("topright", 
         legend=c("Raw counts","Normalized counts"), 
         col=c("red", "blue"), inset = 0.02,box.lty=0,
         cex = 1, lty = 1)
  dev.off()
}



