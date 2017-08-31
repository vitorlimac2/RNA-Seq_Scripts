## read the file

library(ggplot2)
source("/home/vitor/git/RNA-Seq_Scripts/multiplot.R")

### ProR1 vs ProR2

r1_r2 <- read.table("ProR1_ProR2.gene.count.total_nt.mean_nt.nt_by_length")
r1_r2$V2 <- log10(r1_r2$V2*1000000/sum(r1_r2$V2)+1)
r1_r2$V6 <- log10(r1_r2$V6*1000000/sum(r1_r2$V6)+1)



r1_r2$LogFC <- r1_r2$V6 - r1_r2$V2

p1 <- ggplot(r1_r2, aes(x = r1_r2$V2, y = r1_r2$V6 )) + 
  geom_point(aes(colour=LogFC)) + 
  geom_smooth(method='lm',formula=y~x) +
  labs(x = "ProR1", y = "ProR2")
## ProR1 vs ProC2

r1_c2 <- read.table("ProR1_ProC2.gene.count.total_nt.mean_nt.nt_by_length")
cor(r1_c2$V6,r1_c2$V2)

x <- r1_c2[r1_c2$V2>0 & r1_c2$V6 > 0,]

r1_c2 <- x
cor(r1_c2$V6,r1_c2$V2)

r1_c2$V2 <- r1_c2$V2*1000000/sum(r1_c2$V2)
r1_c2$V6 <- r1_c2$V6*1000000/sum(r1_c2$V6)


r1_c2$LogFC <- r1_c2$V2/r1_c2$V6

r1_c2$V2 <- log10(r1_c2$V2*1000000/sum(r1_c2$V2))
r1_c2$V6 <- log10(r1_c2$V6*1000000/sum(r1_c2$V6))

cor(r1_c2$V6,r1_c2$V2)

ggplot(r1_c2, aes(x = r1_c2$V2, 
                        y = r1_c2$V6))+
  geom_point(aes(size=r1_c2$LogFC)) + 
  geom_smooth(method='lm',formula=y~x) + 
  labs(x = "ProR1", y = "ProC2")



r1_c2$V3 <- log10(r1_c2$V3*1000000/sum(r1_c2$V3)+1)
r1_c2$V7 <- log10(r1_c2$V7*1000000/sum(as.numeric(r1_c2$V7))+1)

cor(r1_c2$V7,r1_c2$V3)

r1_c2$nt_counts_r1 <- log10(r1_c2$V2*1000000/sum(r1_c2$V3)+1)
r1_c2$nt_counts_c2 <- log10(r1_c2$V6*1000000/sum(as.numeric(r1_c2$V7))+1)

r1_c2$LogFC <- r1_c2$V6 - r1_c2$V2

r1_c2$LogFC_nt <- r1_c2$V7 - r1_c2$V3

r1_c2$LogFC_new <- r1_c2$nt_counts_c2 - r1_c2$nt_counts_r1

p2 <- ggplot(r1_c2, aes(x = r1_c2$V2, 
                  y = r1_c2$V6))+
  geom_point(aes(colour=LogFC)) + 
  geom_smooth(method='lm',formula=y~x) + 
  labs(x = "ProR1", y = "ProC2")

p2_1 <- ggplot(r1_c2, aes(x = r1_c2$V3, 
                          y = r1_c2$V7))+
  geom_point(aes(colour=LogFC_nt)) + 
  geom_smooth(method='lm',formula=y~x) + 
  labs(x = "ProR1", y = "ProC2")

ggplot(r1_c2, aes(x = r1_c2$nt_counts_r1, 
                          y = r1_c2$nt_counts_c2))+
  geom_point(aes(colour=LogFC_new)) + 
  geom_smooth(method='lm',formula=y~x) + 
  labs(x = "ProR1", y = "ProC2")




multiplot(p2, p2_1, p2_2, cols=2)

## ProR1 vs ProC3

r1_c3 <- read.table("ProR1_ProC3.gene.count.total_nt.mean_nt.nt_by_length")
r1_c3$V2 <- log10(r1_c3$V2*1000000/sum(r1_c3$V2)+1)
r1_c3$V6 <- log10(r1_c3$V6*1000000/sum(r1_c3$V6)+1)

r1_c3$LogFC <- r1_c3$V6 - r1_c3$V2

p3 <- ggplot(r1_c3, aes(x = r1_c3$V2, 
                  y = r1_c3$V6 ))+
  geom_point(aes(colour=LogFC)) + 
  geom_smooth(method='lm',formula=y~x) +
  labs(x = "ProR1", y = "ProC3")



## ProR2 vs ProC2

r2_c2 <- read.table("ProR2_ProC2.gene.count.total_nt.mean_nt.nt_by_length")
r2_c2$V2 <- log10(r2_c2$V2*1000000/sum(r2_c2$V2)+1)
r2_c2$V6 <- log10(r2_c2$V6*1000000/sum(r2_c2$V6)+1)
r2_c2$LogFC <- r2_c2$V6 - r2_c2$V2
p4 <- ggplot(r2_c2, aes(x = r2_c2$V2, 
                  y = r2_c2$V6 ))+
  geom_point(aes(colour=LogFC)) + 
  geom_smooth(method='lm',formula=y~x) +
  labs(x = "ProR2", y = "ProC2")




## ProR2 vs ProC3

r2_c3 <- read.table("ProR2_ProC3.gene.count.total_nt.mean_nt.nt_by_length")
r2_c3$V2 <- log10(r2_c3$V2*1000000/sum(r2_c3$V2)+1)
r2_c3$V6 <- log10(r2_c3$V6*1000000/sum(r2_c3$V6)+1)
r2_c3$LogFC <- r2_c3$V6 - r2_c3$V2
p5 <- ggplot(r2_c3, aes(x = r2_c3$V2, y = r2_c3$V6 ))+
  geom_point(aes(colour=LogFC)) + 
  geom_smooth(method='lm',formula=y~x)+
  labs(x = "ProR2", y = "ProC3")


source("/home/vitor/git/RNA-Seq_Scripts/multiplot.R")

multiplot(p1, p2, p3, p4, p5, cols=2)

###################################################################

r1_c2 <- r1_c2[r1_c2$V2 > 0 & r1_c2$V6 > 0,]

r1r2 <- r1_r2[r1_r2$V2 > 0 & r1_r2$V6 > 0,]

r1c2 <- cbind(r1_c2$V2, r1_c2$V6)

r1c2_ordered <- r1c2[order(r1c2[,2], decreasing =T),]

r1c2_ordered[,1] <- r1c2_ordered[,1]*1000000/sum(r1c2_ordered[,1])
r1c2_ordered[,2] <- r1c2_ordered[,2]*1000000/sum(r1c2_ordered[,2])

cor(r1c2_ordered)

a = 0;
b = nrow(r1c2_ordered); 

plot(r1c2_ordered[a:b,1],r1c2_ordered[a:b,2], log="xy")

