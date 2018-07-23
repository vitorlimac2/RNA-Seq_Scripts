##### Differences of domains
setwd("/home/vitor/final_results/gencode/")
diff_gencode <- read.table("gencode_diff_clusters")$V1
diff_fly <- read.table("dm3_diff_clusters")$V1
diff_worm <- read.table("ce6_diff_clusters")$V1

# my_breaks = pretty(range(c(diff_fly, diff_gencode, diff_worm)), n=10)
# 
# h_gencode <- hist(diff_gencode, breaks=my_breaks, plot=FALSE)
# h_fly <- hist(diff_fly, breaks=my_breaks, plot=FALSE)
# h_worm <- hist(diff_worm, breaks=my_breaks, plot=FALSE)
# 
# y_gencode <- h_gencode$counts*100/length(diff_gencode)
# y_fly <- h_fly$counts*100/length(diff_fly)
# y_worm <- h_worm$counts*100/length(diff_worm)
# 
# max_y <- max(y_gencode, y_fly, y_worm)
# 
# plot(h_gencode$mids, y_gencode, type="l", xlab = "Bitscore/MaxBitscore", ylab = "% differences", col="black", ylim=c(0,max_y)) 

my_breaks = pretty(range(c(diff_gencode)), n=10)

h_gencode <- hist(diff_gencode, breaks=my_breaks, plot=FALSE)

y_gencode <- h_gencode$counts*100/length(diff_gencode)

max_y <- max(y_gencode)

########################

barplot(col="white", y_gencode, axes = F, space = 0, width = 0.2, ylim = c(0,30))
axis(2, at=seq(0,30,by=5))
labels = c("0","<0.1","<0.2","<0.3","<0.4","<0.5","<0.6","<0.7","<0.8","<0.9","<=1")
axis(1, at=seq(-0.1, 2, by=0.20), labels = labels)
title(ylab = "Domain cluster [%]", xlab = "Bitscore/MaxScore of each cluster")



x = h_gencode$mids
y = y_gencode
smoothingSpline = smooth.spline(x, y, spar=0.35)
plot(x,y, col="darkgray", lwd=3, ylab="Domain Clusters [%]", xlab="Bitscores/MaxScore of each cluster", xaxt='n')
lines(smoothingSpline)
axis(1, at=seq(0,1,by=0.1))


# x = h_fly$mids
# y = y_fly
# smoothingSpline = smooth.spline(x, y, spar=0.35)
# lines(smoothingSpline, col="black", lwd=3,)
# 
# x = h_worm$mids
# y = y_worm
# smoothingSpline = smooth.spline(x, y, spar=0.35)
# lines(smoothingSpline, col="darkgray", lwd=3,)
# 

################################################

gencode_L_diff <- read.table("gencode_hmmLength_clusterDiff")

dlengths10 <- gencode_L_diff[gencode_L_diff$V2 < 0.1,]$V1
dlengths20 <- gencode_L_diff[gencode_L_diff$V2 >= 0.1 & gencode_L_diff$V2 < 0.2,]$V1
dlengths30 <- gencode_L_diff[gencode_L_diff$V2 >= 0.2 & gencode_L_diff$V2 < 0.3,]$V1
dlengths40 <- gencode_L_diff[gencode_L_diff$V2 >= 0.3 & gencode_L_diff$V2 < 0.4,]$V1
dlengths50 <- gencode_L_diff[gencode_L_diff$V2 >= 0.4 & gencode_L_diff$V2 < 0.5,]$V1
dlengths60 <- gencode_L_diff[gencode_L_diff$V2 >= 0.5 & gencode_L_diff$V2 < 0.6,]$V1
dlengths70 <- gencode_L_diff[gencode_L_diff$V2 >= 0.6 & gencode_L_diff$V2 < 0.7,]$V1
dlengths80 <- gencode_L_diff[gencode_L_diff$V2 >= 0.7 & gencode_L_diff$V2 < 0.8,]$V1
dlengths90 <- gencode_L_diff[gencode_L_diff$V2 >= 0.8 & gencode_L_diff$V2 < 0.9,]$V1
dlengths100 <- gencode_L_diff[gencode_L_diff$V2 >= 0.9 & gencode_L_diff$V2 <= 1,]$V1

dd<- lapply(c("dlengths10", "dlengths20", "dlengths30", "dlengths40","dlengths50","dlengths60",
           "dlengths70","dlengths80","dlengths90","dlengths100"),get,envir=environment())

boxplot(dd, outline = F, ylab = "Domain lengths (AA)", axes = F)
axis(2, at=seq(0,900, by=100))
labels = c("0","<0.1","<0.2","<0.3","<0.4","<0.5","<0.6","<0.7","<0.8","<0.9","<=1")
axis(1, at=0:10, labels = labels)

########################################################################

setwd("/home/vitor/final_results/gencode")
gencode_L_diff <- read.table("gencode_numPred_clusterDiff")

dlengths10 <- gencode_L_diff[gencode_L_diff$V2 < 0.1,]$V1
dlengths20 <- gencode_L_diff[gencode_L_diff$V2 >= 0.1 & gencode_L_diff$V2 < 0.2,]$V1
dlengths30 <- gencode_L_diff[gencode_L_diff$V2 >= 0.2 & gencode_L_diff$V2 < 0.3,]$V1
dlengths40 <- gencode_L_diff[gencode_L_diff$V2 >= 0.3 & gencode_L_diff$V2 < 0.4,]$V1
dlengths50 <- gencode_L_diff[gencode_L_diff$V2 >= 0.4 & gencode_L_diff$V2 < 0.5,]$V1
dlengths60 <- gencode_L_diff[gencode_L_diff$V2 >= 0.5 & gencode_L_diff$V2 < 0.6,]$V1
dlengths70 <- gencode_L_diff[gencode_L_diff$V2 >= 0.6 & gencode_L_diff$V2 < 0.7,]$V1
dlengths80 <- gencode_L_diff[gencode_L_diff$V2 >= 0.7 & gencode_L_diff$V2 < 0.8,]$V1
dlengths90 <- gencode_L_diff[gencode_L_diff$V2 >= 0.8 & gencode_L_diff$V2 < 0.9,]$V1
dlengths100 <- gencode_L_diff[gencode_L_diff$V2 >= 0.9 & gencode_L_diff$V2 <= 1,]$V1

dd<- lapply(c("dlengths10", "dlengths20", "dlengths30", "dlengths40","dlengths50","dlengths60",
              "dlengths70","dlengths80","dlengths90","dlengths100"),get,envir=environment())

boxplot(dd, col = "gray", outline = F, ylab = "Domain predictions", axes = F)
axis(2, at=seq(0,350, by=50))
labels = c("0","<0.1","<0.2","<0.3","<0.4","<0.5","<0.6","<0.7","<0.8","<0.9","<=1")
axis(1, at=0:10, labels = labels)

