r1 <- read.table("ProR_1.fastq.trimmed.fq.read.lengths")$V1
r2 <- read.table("ProR_2.fastq.trimmed.fq.read.lengths")$V1
c1 <- read.table("ProC_1.fastq.trimmed.fq.read.lengths")$V1
c2 <- read.table("ProC_2.fastq.trimmed.fq.read.lengths")$V1
c3 <- read.table("ProC_3.fastq.trimmed.fq.read.lengths")$V1
c4 <- read.table("ProC_4.fastq.trimmed.fq.read.lengths")$V1
c5 <- read.table("ProC_5.fastq.trimmed.fq.read.lengths")$V1
c6 <- read.table("ProC_6.fastq.trimmed.fq.read.lengths")$V1
c7 <- read.table("ProC_7.fastq.trimmed.fq.read.lengths")$V1
c8 <- read.table("ProC_8.fastq.trimmed.fq.read.lengths")$V1
c9 <- read.table("ProC_9.fastq.trimmed.fq.read.lengths")$V1

reps <- c("R1","R2","C1","C2","C3","C4","C5","C6","C7","C8","C9")

#boxplot(r1,r2,c1,c2,c3,c4,c5,c6,c7,c8,c9, names=reps)


my_breaks <- pretty(range(c(r1,r2, c1,c2,c3,
                            c4,c5,c6,
                            c7,c8,c9)), n=15)

h_r1 <- hist(r1, breaks=my_breaks, plot=FALSE)
h_r2 <- hist(r2, breaks=my_breaks, plot=FALSE)
h_c1 <- hist(c1, breaks=my_breaks, plot=FALSE)
h_c2 <- hist(c2, breaks=my_breaks, plot=FALSE)
h_c3 <- hist(c3, breaks=my_breaks, plot=FALSE)
h_c4 <- hist(c4, breaks=my_breaks, plot=FALSE)
h_c5 <- hist(c5, breaks=my_breaks, plot=FALSE)
h_c6 <- hist(c6, breaks=my_breaks, plot=FALSE)
h_c7 <- hist(c7, breaks=my_breaks, plot=FALSE)
h_c8 <- hist(c8, breaks=my_breaks, plot=FALSE)
h_c9 <- hist(c9, breaks=my_breaks, plot=FALSE)


max_y <- max(h_r1$counts*100/sum(length(r1)),
             h_r2$counts*100/sum(length(r2)),
             h_c1$counts*100/sum(length(c1)),
             h_c2$counts*100/sum(length(c2)),
             h_c3$counts*100/sum(length(c3)),
             h_c4$counts*100/sum(length(c4)),
             h_c5$counts*100/sum(length(c5)),
             h_c6$counts*100/sum(length(c6)),
             h_c7$counts*100/sum(length(c7)),
             h_c8$counts*100/sum(length(c8)),
             h_c9$counts*100/sum(length(c9)))

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")

lines(h_r1$mids, h_r1$counts*100/sum(length(r1)), type="l", ylab = "%", xlab = "Read length (bp)", col="blue", ylim=c(0,max_y), xlim=c(35,250)) 

lines(h_r2$mids, h_r2$counts*100/sum(length(r2)), col="green")
lines(h_c1$mids, h_c1$counts*100/sum(length(c1)), col="purple")
lines(h_c2$mids, h_c2$counts*100/sum(length(c2)), col="black")
lines(h_c3$mids, h_c3$counts*100/sum(length(c3)), col="red")
lines(h_c4$mids, h_c4$counts*100/sum(length(c4)), col="yellow")
lines(h_c5$mids, h_c5$counts*100/sum(length(c5)), col= "orange")
lines(h_c6$mids, h_c6$counts*100/sum(length(c6)), col="brown")
lines(h_c7$mids, h_c7$counts*100/sum(length(c7)), col="darkgoldenrod4")
lines(h_c8$mids, h_c8$counts*100/sum(length(c8)), col="darkgreen")
lines(h_c9$mids, h_c9$counts*100/sum(length(c9)), col = "pink")

par(bg = 'white')

legend(230,35, cex=0.8, reps, fill=c("blue",
                                     "green",
                                     "purple",
                                     "black",
                                     "red",
                                     "yellow",
                                     "orange",
                                     "brown",
                                     "darkgoldenrod4",
                                     "darkgreen",
                                     "pink"), horiz=FALSE)


