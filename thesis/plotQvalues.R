workdir <- "/home/vitor/Downloads"
setwd(workdir)
pror1 <- read.table("ProR1.GeneDistFrag.ks.output.Qvalues", skip = 1, row.names = 1)
colnames(pror1) <- c("less_distance", "less", "twosided_distance","twosided","greater_distance","greater")
head(pror1)

my_breaks <- pretty(range(c(pror1$less, pror1$twosided, pror1$greater)), n=10)

h_l <- hist(pror1$less, breaks=my_breaks, plot=FALSE)
h_t <- hist(pror1$twosided, breaks=my_breaks, plot=FALSE)
h_g <- hist(pror1$greater, breaks=my_breaks, plot=FALSE)

max_y <- max(h_l$counts*100/length(pror1$less), h_t$counts*100/length(pror1$twosided),
             h_g$counts*100/length(pror1$greater))

plot(h_l$mids, h_l$counts*100/length(pror1$less), type="l", ylab = "%", xlab = "Q-values", col="blue", ylim=c(0,max_y), xlim=c(0,1))

lines(h_t$mids, h_t$counts*100/length(pror1$twosided), col="red")

lines(h_g$mids, h_g$counts*100/length(pror1$twosided), col="black")

boxplot(pror1[pror1$less < 0.01,]$less_distance, outline = F)
boxplot(pror1[pror1$twosided < 0.01,]$less_distance, outline = F)
boxplot(pror1[pror1$greater < 0.01,]$less_distance, outline = F)

summary(pror1[pror1$less < 0.01 & pror1$less_distance > 0.5,]$less_distance)
