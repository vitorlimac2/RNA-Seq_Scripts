workdir <- "/home/vitor/Proj_ProC_R/mappings_star/"

replicates <- c("ProR1",
                "ProR2")

r1_less <- log(read.table(paste(workdir,"ProR1.ksTest.Less.Significant.MedianTx.GC_lengths",sep=""))$V1)
r1_greater <- log(read.table(paste(workdir,"ProR1.ksTest.Greater.Significant.MedianTx.GC_lengths",sep=""))$V1)
r2_less <- log(read.table(paste(workdir,"ProR2.ksTest.Less.Significant.MedianTx.GC_lengths",sep=""))$V1)
r2_greater <- log(read.table(paste(workdir,"ProR2.ksTest.Greater.Significant.MedianTx.GC_lengths",sep=""))$V1)


my_breaks <- pretty(range(c(r1_less, r1_greater,r2_less, r2_greater)), n=20)




h_r1_less <- hist(r1_less, breaks=my_breaks, plot=FALSE)
h_r1_greater <- hist(r1_greater, breaks=my_breaks, plot=FALSE)
h_r2_less <- hist(r2_less, breaks=my_breaks, plot=FALSE)
h_r2_greater <- hist(r1_greater, breaks=my_breaks, plot=FALSE)

max_y <- max(h_r1_less$counts*100/length(r1_less),
             h_r1_greater$counts*100/length(r1_greater),
             h_r2_less$counts*100/length(r2_less),
             h_r2_greater$counts*100/length(r2_greater))


## LESS group - Genes with more long reads
## GREATER group - gene with more short reads


plot(main="ProR1", lty=2,h_r1_less$mids, h_r1_less$counts*100/length(r1_less), type="l", ylab = "%", xlab = "Median Transcript Length (log)", col="black", ylim=c(0,max_y)) 
lines(h_r1_greater$mids, h_r1_greater$counts*100/length(r1_greater), col="grey")

legend("topright",lwd=2, cex=0.8, bty = "n", legend=c("Long Reads", "Short Reads"), col=c("black","grey"),
       lty = c(2,1),
       horiz=FALSE)

plot(main="ProR2", lty=2,h_r2_less$mids, h_r2_less$counts*100/length(r2_less), type="l", ylab = "%", xlab = "Median Transcript Length (log)", col="black", ylim=c(0,max_y)) 
lines(h_r2_greater$mids, h_r2_greater$counts*100/length(r2_greater), col="grey")

legend("topright",lwd=2, cex=0.8, bty = "n", legend=c("Long Reads", "Short Reads"), col=c("black","grey"),
       lty = c(2,1),
       horiz=FALSE)

#################################################################################################
#################################################################################################

r1_less <- read.table(paste(workdir,"ProR1.ksTest.Less.Significant.MedianTx.GC_lengths",sep=""),fill = !blank.lines.skip)$V2
r1_greater <- read.table(paste(workdir,"ProR1.ksTest.Greater.Significant.MedianTx.GC_lengths",sep=""),fill = !blank.lines.skip)$V2
r2_less <- read.table(paste(workdir,"ProR2.ksTest.Less.Significant.MedianTx.GC_lengths",sep=""),fill = !blank.lines.skip)$V2
r2_greater <- read.table(paste(workdir,"ProR2.ksTest.Greater.Significant.MedianTx.GC_lengths",sep=""),fill = !blank.lines.skip)$V2


my_breaks <- pretty(range(c(r1_less, r1_greater,r2_less, r2_greater)), n=20)




h_r1_less <- hist(r1_less, breaks=my_breaks, plot=FALSE)
h_r1_greater <- hist(r1_greater, breaks=my_breaks, plot=FALSE)
h_r2_less <- hist(r2_less, breaks=my_breaks, plot=FALSE)
h_r2_greater <- hist(r1_greater, breaks=my_breaks, plot=FALSE)

max_y <- max(h_r1_less$counts*100/length(r1_less),
             h_r1_greater$counts*100/length(r1_greater),
             h_r2_less$counts*100/length(r2_less),
             h_r2_greater$counts*100/length(r2_greater))


## LESS group - Genes with more long reads
## GREATER group - gene with more short reads


plot(main="ProR1", lty=2,h_r1_less$mids, h_r1_less$counts*100/length(r1_less), type="l", ylab = "% Genes (median transcript length)", xlab = "% GC", col="black", xlim=c(20,80), ylim=c(0,max_y)) 
lines(h_r1_greater$mids, h_r1_greater$counts*100/length(r1_greater), col="grey")

legend("topright",lwd=2, cex=0.8, bty = "n", legend=c("Long Reads", "Short Reads"), col=c("black","grey"),
       lty = c(2,1),
       horiz=FALSE)

plot(main="ProR2", lty=2,h_r2_less$mids, h_r2_less$counts*100/length(r2_less), type="l", ylab = "% Genes (median transcript length)", xlab = "% GC", col="black", ylim=c(0,max_y), xlim=c(20,80)) 
lines(h_r2_greater$mids, h_r2_greater$counts*100/length(r2_greater), col="grey")

legend("topright",lwd=2, cex=0.8, bty = "n", legend=c("Long Reads", "Short Reads"), col=c("black","grey"),
       lty = c(2,1),
       horiz=FALSE)