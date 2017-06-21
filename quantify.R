####################### MULTIMAPPINGS ############################################

proR1.table <-read.table("ProR1.sorted.by.median.length.counts") 
proC5.table <- read.table("ProC5.sorted.by.median.length.counts")

proR1_50.mm <- read.table("../ProR_1_50_sorted.bam.multimaps")
proR1_100.mm <- read.table("../ProR_1_100_sorted.bam.multimaps")
proR1_150.mm <- read.table("../ProR_1_150_sorted.bam.multimaps")
proR1_200.mm <- read.table("../ProR_1_200_sorted.bam.multimaps")
proR1_250.mm <- read.table("../ProR_1_250_sorted.bam.multimaps")
proR1_300.mm <- read.table("../ProR_1_300_sorted.bam.multimaps")

proR1.mm <- c(proR1_50.mm$V1, proR1_100.mm$V1, proR1_150.mm$V1, proR1_200.mm$V1, proR1_250.mm$V1, proR1_300.mm$V1)

hist_proR1_mm <- hist(proR1.mm, plot = F, breaks = c(0,1,2,3,4))

hist_proR1_mm$counts <- log10(hist_proR1_mm$counts)

plot( hist_proR1_mm$breaks[0:4], hist_proR1_mm$counts, 
      ylab='log10(Frequency)', type="o", col="blue", xaxt='n', main="Multimappings",xlab="# of Multimap",
      ylim=c(0,7))

axis(1, at=c(0,1,2,3), labels=c(1,2,3,4))

proC5_50.mm <- read.table("../ProC_5_50_sorted.bam.multimaps")
proC5_100.mm <- read.table("../ProC_5_100_sorted.bam.multimaps")
proC5_150.mm <- read.table("../ProC_5_150_sorted.bam.multimaps")
proC5_200.mm <- read.table("../ProC_5_200_sorted.bam.multimaps")
proC5_250.mm <- read.table("../ProC_5_250_sorted.bam.multimaps")
proC5_300.mm <- read.table("../ProC_5_300_sorted.bam.multimaps")

proC5.mm <- c(proC5_50.mm$V1, proC5_100.mm$V1, proC5_150.mm$V1, proC5_200.mm$V1, 
              proC5_250.mm$V1, proC5_300.mm$V1)

hist_proC5_mm <- hist(proC5.mm, plot = F, breaks = c(0,1,2,3,4))

hist_proC5_mm$counts <- log10(hist_proC5_mm$counts)

lines(hist_proC5_mm$breaks[0:4], hist_proC5_mm$counts, col="red", type = "o")

legend("topright", inset=.05, c("ProR1","ProC5"), fill=c("blue","red"), horiz=TRUE)

total_proR1 <- sum(hist_proR1_mm$counts)
total_proC5 <- sum(hist_proC5_mm$counts)

plot( hist_proR1_mm$breaks[0:4], hist_proR1_mm$counts*100/total_proR1, 
      ylab='%', type="o", col="blue", xaxt='n', main="Multimappings",xlab="# of Multimap",
      ylim=c(0,45))

axis(1, at=c(0,1,2,3), labels=c(1,2,3,4))

lines(hist_proC5_mm$breaks[0:4], hist_proC5_mm$counts*100/total_proC5, col="red", type = "o")

proR1_mm_length <- c(sum(proR1_50.mm$V1), sum(proR1_100.mm$V1), sum(proR1_150.mm$V1),
                     sum(proR1_200.mm$V1),sum(proR1_250.mm$V1), sum(proR1_300.mm$V1))*100/sum(proR1.mm)


proC5_mm_length <- c(sum(proC5_50.mm$V1), sum(proC5_100.mm$V1), sum(proC5_150.mm$V1),
                     sum(proC5_200.mm$V1),sum(proC5_250.mm$V1), sum(proC5_300.mm$V1))*100/sum(proC5.mm)

plot( proR1_mm_length, 
      ylab='%', type="o", col="blue", xaxt='n', main="Multimappings",xlab="Read length (bp)",
      ylim=c(0,80))

axis(1, at=c(1, 2, 3 , 4, 5, 6), labels=c("[20,50[", "[50,100[",
                                          "[100,150[","[150,200[",
                                          "[200,250[",">=250"))
lines(c(1:6), proC5_mm_length, col="red", type = "o")

legend("topleft", inset=.05, c("ProR1","ProC5"), fill=c("blue","red"), horiz=TRUE)

########################### GENE MEDIAN LENGTHS VS READ LENGTHS #######################


# proR1.50.genes.read <- read.table("ProR_1_50_geneLength_readLength")
# proR1.100.genes.read <- read.table("ProR_1_100_geneLength_readLength")
# proR1.150.genes.read <- read.table("ProR_1_150_geneLength_readLength")
# proR1.200.genes.read <- read.table("ProR_1_200_geneLength_readLength")
# proR1.250.genes.read <- read.table("ProR_1_250_geneLength_readLength")
# proR1.300.genes.read <- read.table("ProR_1_300_geneLength_readLength")
# 
# 
# 
# proC5.50.genes.read <- read.table("ProC_5_50_geneLength_readLength")
# proC5.100.genes.read <- read.table("ProC_5_100_geneLength_readLength")
# proC5.150.genes.read <- read.table("ProC_5_150_geneLength_readLength")
# proC5.200.genes.read <- read.table("ProC_5_200_geneLength_readLength")
# proC5.250.genes.read <- read.table("ProC_5_250_geneLength_readLength")
# proC5.300.genes.read <- read.table("ProC_5_300_geneLength_readLength")

proR1.genes.read <-read.table("ProR_1_geneLength_readLength")

proC5.genes.read <-read.table("ProC_5_geneLength_readLength")

summary(proR1.genes.read)
summary(proC5.genes.read)
quantile(x = proR1.genes.read$V1, probs = seq(0,1,0.1))
quantile(x = proC5.genes.read$V1, probs = seq(0,1,0.1))

##################### 10% shorter and longer
### ProR1
## 10% shorter <582
## 10% longer >2582

### 644 773 962 1685

### ProC5
## 10% shorter < 582
## 10% longer >2671

### 644 771 956 1706
ProR1_set20 <- proR1.genes.read[proR1.genes.read$V1<644,]$V2
ProR1_set40 <- proR1.genes.read[proR1.genes.read$V1>=644 & proR1.genes.read$V1<773,]$V2
ProR1_set60 <- proR1.genes.read[proR1.genes.read$V1>=733 & proR1.genes.read$V1<962,]$V2
ProR1_set80 <- proR1.genes.read[proR1.genes.read$V1>=962 & proR1.genes.read$V1<1685,]$V2
ProR1_set100 <- proR1.genes.read[proR1.genes.read$V1>=1685,]$V2


ProC5_set20 <- proC5.genes.read[proC5.genes.read$V1<644,]$V2
ProC5_set40 <- proC5.genes.read[proC5.genes.read$V1>=644 & proC5.genes.read$V1<771,]$V2
ProC5_set60 <- proC5.genes.read[proC5.genes.read$V1>=771 & proC5.genes.read$V1<956,]$V2
ProC5_set80 <- proC5.genes.read[proC5.genes.read$V1>=956 & proC5.genes.read$V1<1706,]$V2
ProC5_set100 <- proC5.genes.read[proC5.genes.read$V1>=1706,]$V2

boxplot(ylab="Mapped Read Length (bp)", 
        ProR1_set20,
        ProC5_set20,
        ProR1_set40,
        ProC5_set40,
        ProR1_set60,
        ProC5_set60,
        ProR1_set80,
        ProC5_set80,
        ProR1_set100,
        ProC5_set100,
        col = c("white","gray"),
        outline = F, xaxt="n",
        xlab="Median Transcript Length (Quantile probs)")

abline(v=c(2.5, 4.5, 6.5, 8.5))

axis(1,at=seq(1.5,9.5,by=2), labels=c("0-20%","20-40%","40-60%","60-80%","80-100%"))

legend("topleft", inset=.002, c("ProR1","ProC5"), fill=c("white","gray"), cex=0.8)


ProR1_set20 <- proR1.genes.read[proR1.genes.read$V1<644,]$V2
ProR1_set40 <- proR1.genes.read[proR1.genes.read$V1>=644 & proR1.genes.read$V1<773,]$V2
ProR1_set60 <- proR1.genes.read[proR1.genes.read$V1>=733 & proR1.genes.read$V1<962,]$V2
ProR1_set80 <- proR1.genes.read[proR1.genes.read$V1>=962 & proR1.genes.read$V1<1685,]$V2
ProR1_set100 <- proR1.genes.read[proR1.genes.read$V1>=1685,]$V2


ProC5_set20 <- proC5.genes.read[proC5.genes.read$V1<644,]$V2
ProC5_set40 <- proC5.genes.read[proC5.genes.read$V1>=644 & proC5.genes.read$V1<771,]$V2
ProC5_set60 <- proC5.genes.read[proC5.genes.read$V1>=771 & proC5.genes.read$V1<956,]$V2
ProC5_set80 <- proC5.genes.read[proC5.genes.read$V1>=956 & proC5.genes.read$V1<1706,]$V2
ProC5_set100 <- proC5.genes.read[proC5.genes.read$V1>=1706,]$V2

###############################################################################
##### Nucleotide Coverage vs. Read Coverage

nut_ProR1_set20 <- sum(proR1.genes.read[proR1.genes.read$V1<644,]$V2)
nut_ProR1_set40 <- sum(proR1.genes.read[proR1.genes.read$V1>=644 & proR1.genes.read$V1<773,]$V2)
nut_ProR1_set60 <- sum(proR1.genes.read[proR1.genes.read$V1>=733 & proR1.genes.read$V1<962,]$V2)
nut_ProR1_set80 <- sum(proR1.genes.read[proR1.genes.read$V1>=962 & proR1.genes.read$V1<1685,]$V2)
nut_ProR1_set100 <- sum(proR1.genes.read[proR1.genes.read$V1>=1685,]$V2)


nut_ProR1_vector <- c(nut_ProR1_set20, nut_ProR1_set40, nut_ProR1_set60,
                      nut_ProR1_set80,nut_ProR1_set100)

read_ProR1_vector <- c(length(ProR1_set20), length(ProR1_set40), length(ProR1_set60),
                       length(ProR1_set80), length(ProR1_set100))


###############################################################################
###### Gene counts vs Length

proR1_count_reads <- read.table("ProR_1_count_readLength")

proC5_count_reads <- read.table("ProC_5_count_readLength")


quantile(x = proR1_count_reads$V1, probs = seq(0,1,0.1))

quantile(x = proC5_count_reads$V1, probs = seq(0,1,0.1))

### comprimento dos 20% menos/mais expressos

proR1_count_set20 <- proR1_count_reads[proR1_count_reads$V1 < 129,]$V2
proC5_count_set20 <- proC5_count_reads[proC5_count_reads$V1 < 172.4,]$V2

proR1_count_set100 <- proR1_count_reads[proR1_count_reads$V1 > 3230,]$V2
proC5_count_set100 <- proC5_count_reads[proC5_count_reads$V1 > 4703,]$V2

boxplot(proR1_count_set20,
        proC5_count_set20,
        proR1_count_set100,
        proC5_count_set100, outline = F)

######################################################
##### Faixa de transcritos

ProR1_set200 <- proR1.genes.read[proR1.genes.read$V1<200,]$V2
ProR1_set1000 <- proR1.genes.read[proR1.genes.read$V1>=900 & proR1.genes.read$V1<1000,]$V2
ProR1_set2000 <- proR1.genes.read[proR1.genes.read$V1>=10000,]$V2

ProC5_set200 <- proC5.genes.read[proC5.genes.read$V1<200,]$V2
ProC5_set1000 <- proC5.genes.read[proC5.genes.read$V1>=900 & proC5.genes.read$V1<1000,]$V2
ProC5_set2000 <- proC5.genes.read[proC5.genes.read$V1>=10000,]$V2



d1 <- ProR1_set200
d3 <- ProC5_set200 

breaks <- pretty(range(c(d1, d3)), n=20)

D1 <- hist(d1, breaks=breaks, plot=FALSE)$counts*100/sum(d1)
D3 <- hist(d3, breaks=breaks, plot=FALSE)$counts*100/sum(d3)


#dat <- rbind(D1, D2)
dat <- rbind(D1, D3)
colnames(dat) <- paste(breaks[-length(breaks)], breaks[-1], sep="-")

barplot(dat, beside=TRUE, space=c(0, 0.1), las=2, main="Transcript length < 200\n
        ProR1 vs ProC5")

legend("topright", inset=.05, c("ProR1","ProC5"), fill=c("dimgrey","azure"), horiz=F)


d2 <- ProR1_set2000
breaks <- pretty(range(c(d1, d2)), n=20)
D2 <- hist(d2, breaks=breaks, plot=FALSE)$counts*100/sum(d2)
D1 <- hist(d1, breaks=breaks, plot=FALSE)$counts*100/sum(d1)
mydat <- rbind(D1, D2)
colnames(mydat) <- paste(breaks[-length(breaks)], breaks[-1], sep="-")

barplot(mydat, beside=TRUE, space=c(0, 0.1), las=2, main="ProR1")

legend("topright", inset=.05, c("length < 200","length > 10,000"), fill=c("dimgrey","azure"), horiz=F)

########################################################

ProC5_set200 <- proC5.genes.read[proC5.genes.read$V1<200,]$V2
ProC5_set2000 <- proC5.genes.read[proC5.genes.read$V1>=10000,]$V2

d1 <- ProC5_set200
d2 <- ProC5_set2000

breaks <- pretty(range(c(d1, d2)), n=20)

D1 <- hist(d1, breaks=breaks, plot=FALSE)$counts
D2 <- hist(d2, breaks=breaks, plot=FALSE)$counts


#dat <- rbind(D1, D2)
dat <- rbind(D1, D2)
colnames(dat) <- paste(breaks[-length(breaks)], breaks[-1], sep="-")

barplot(dat, beside=TRUE, space=c(0, 0.1), las=2, main="ProR1")

legend("topright", title="Transcript length", inset=.05, c("< 200","> 10,000"), fill=c("dimgrey","azure"), horiz=F)

