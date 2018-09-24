setwd("/home/vitor/Proj_ProC_R/analysis/")

rep1 <- read.table("ProR1.MappedReadInfo.Filled", sep = "\t")[,2:3]

rep2 <- read.table("ProC2.MappedReadInfo.Filled", sep = "\t")[,2:3]

rep3 <- read.table("ProC8.MappedReadInfo.Filled", sep = "\t")[,2:3]

rep1.unique.mapped <- rep1[rep1[,2] == 1,]
rep2.unique.mapped <- rep2[rep2[,2] == 1,]
rep3.unique.mapped <- rep3[rep3[,2] == 1,]

breaks <- pretty(range(c(rep1[rep[,3] == 1,1], 
                         rep2[rep[,3] == 1,1], 
                         rep3[rep[,3] == 1,1])), n = 5)


h1 <- hist(rep1[rep1[,2] == 1,1], breaks = breaks, plot=FALSE)
h2 <- hist(rep2[rep2[,2] == 1,1], breaks = breaks, plot=FALSE)
h3 <- hist(rep3[rep3[,2] == 1,1], breaks = breaks, plot=FALSE)

dat <- rbind(h1$counts*100/nrow(rep1),h2$counts*100/nrow(rep2),h3$counts*100/nrow(rep3))

colnames(dat) <- paste(breaks[-length(breaks)], breaks[-1], sep="-")

barplot(t(x), beside = TRUE, las=1, space=c(0, 0.2), 
        ylab = "Unique mapped reads %",
        xlab="Read length (bp)")

legend("topright", 
       legend = c("ProR1", "ProC2", "ProC8"), 
       fill = c("darkgray","gray", "lightgray"),
       cex=0.9,
       border = F,
       bty = 'n')

########## GGPLOT ########################

x <- t(dat)[1:6,]
x <- as.data.frame(x)
colnames(x) <- c("ProR1", "ProC2", "ProC8")

x[,1:3] <- round(x[,1:3], digits = 2)

x$ranges <- rownames(x)
x$ranges <- factor(x$ranges, levels = rownames(x))

p <- x %>% gather(Amostra, Porcentagem, ProR1, ProC2, ProC8) %>% 
  ggplot(aes(x=ranges, y=Porcentagem, colour=Amostra, group=Amostra)) + geom_point() + geom_line() + theme_bw()

p <- p + labs(x="Comprimento de read (pb)", y="Reads mapeados unicamente (%)")

p + geom_text(aes(label=Porcentagem, hjust=1, vjust=-0.5))
