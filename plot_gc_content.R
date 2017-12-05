workdir <- "/home/vitor/Proj_ProC_R/mappings_star/"

replicates <- c("ProC1",
                "ProC2",
                "ProC3",
                "ProC4",
                "ProC5",
                "ProC6",
                "ProC7",
                "ProC8",
                "ProC9",
                "ProR1",
                "ProR2")

c1 <- read.table(gzfile(paste(workdir,"ProC1_GC_Content.gz",sep="")))$V2
c2 <- read.table(gzfile(paste(workdir,"ProC2_GC_Content.gz",sep="")))$V2
c3 <- read.table(gzfile(paste(workdir,"ProC3_GC_Content.gz",sep="")))$V2
c4 <- read.table(gzfile(paste(workdir,"ProC4_GC_Content.gz",sep="")))$V2
c5 <- read.table(gzfile(paste(workdir,"ProC5_GC_Content.gz",sep="")))$V2
c6 <- read.table(gzfile(paste(workdir,"ProC6_GC_Content.gz",sep="")))$V2
c7 <- read.table(gzfile(paste(workdir,"ProC7_GC_Content.gz",sep="")))$V2
c8 <- read.table(gzfile(paste(workdir,"ProC8_GC_Content.gz",sep="")))$V2
c9 <- read.table(gzfile(paste(workdir,"ProC9_GC_Content.gz",sep="")))$V2
r1 <- read.table(gzfile(paste(workdir,"ProR1_GC_Content.gz",sep="")))$V2
r2 <- read.table(gzfile(paste(workdir,"ProR2_GC_Content.gz",sep="")))$V2

my_breaks <- pretty(range(c(c1,c2,c3,c4,c5,c6,c7,c8,c9,r1,r2)), n=10)

h_c1 <- hist(c1, breaks=my_breaks, plot=FALSE)
h_c2 <- hist(c2, breaks=my_breaks, plot=FALSE)
h_c3 <- hist(c3, breaks=my_breaks, plot=FALSE)
h_c4 <- hist(c4, breaks=my_breaks, plot=FALSE)
h_c5 <- hist(c5, breaks=my_breaks, plot=FALSE)
h_c6 <- hist(c6, breaks=my_breaks, plot=FALSE)
h_c7 <- hist(c7, breaks=my_breaks, plot=FALSE)
h_c8 <- hist(c8, breaks=my_breaks, plot=FALSE)
h_c9 <- hist(c9, breaks=my_breaks, plot=FALSE)
h_r1 <- hist(r1, breaks=my_breaks, plot=FALSE)
h_r2 <- hist(r2, breaks=my_breaks, plot=FALSE)

max_y <- max(h_c1$counts*100/length(c1)
             , h_c2$counts*100/length(c2)
             , h_c3$counts*100/length(c3)
             , h_c4$counts*100/length(c4)
             , h_c5$counts*100/length(c5)
             , h_c6$counts*100/length(c6)
             , h_c7$counts*100/length(c7)
             , h_c8$counts*100/length(c8)
             , h_c9$counts*100/length(c9)
             , h_r1$counts*100/length(r1)
             , h_r2$counts*100/length(r2))

plot(h_c1$mids, h_c1$counts*100/length(c1), type="l", ylab = "% reads", xlab = "% GC", col="purple", ylim=c(0,max_y), xlim=c(0,100)) 

lines(h_c2$mids, h_c2$counts*100/length(c2), col="black")
lines(h_c3$mids, h_c3$counts*100/length(c3), col="red")
lines(h_c4$mids, h_c4$counts*100/length(c4), col="yellow")
lines(h_c5$mids, h_c5$counts*100/length(c5), col= "orange")
lines(h_c6$mids, h_c6$counts*100/length(c6), col="brown",type = "o")
lines(h_c7$mids, h_c7$counts*100/length(c7), col="darkgoldenrod4",type = "o")
lines(h_c8$mids, h_c8$counts*100/length(c8), col="darkgreen",type = "o")
lines(h_c9$mids, h_c9$counts*100/length(c9), col = "pink",type = "o")

lines(h_r1$mids, h_r1$counts*100/length(r1), col="blue",lty=2) 
lines(h_r2$mids, h_r2$counts*100/length(r2), col="green", lty=2)

legend(80,30, lwd=2, cex=0.8, bty = "n", legend=replicates, col=c("purple",
                                     "black",
                                     "red",
                                     "yellow",
                                     "orange",
                                     "brown",
                                     "darkgoldenrod4",
                                     "darkgreen",
                                     "pink",
                                     "blue",
                                     "green"),
        lty = c(1,1,1,1,1,NA,NA,NA,NA,2,2),
       pch = c(NA,NA,NA,NA,NA,"o","o","o","o",NA,NA),
       horiz=FALSE)
