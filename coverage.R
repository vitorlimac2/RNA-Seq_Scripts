setwd("/home/vitor/PRJ.SRP064142/BWA_ION_QUANT/2-quantification-featureCounts/")

#####################################################################

cov_nt <- read.table("ProR1_ProC5.gene.median_length.count.total_nt.mean_nt.nt_by_length",row.names = 1)

head(cov_nt)

colnames(cov_nt) <- c("length",
                      "R1_count", 
                      "R1_nt_total",
                      "R1_nt_mean",
                      "R1_nt_by_length",
                      "C5_count", 
                      "C5_nt_total",
                      "C5_nt_mean",
                      "C5_nt_by_length")
detach(cov_nt)
attach(cov_nt)

p <- ggplot(cov_nt, aes(x = log10(R1_count+1), 
                                y = log10(C5_count+1)),
                        color = R1_nt_mean) +
  geom_point(alpha=0.2) +
  stat_smooth(method="lm", se=T, col="red")
  

p + annotate("text", label = "R² = 0.8987; p-value < 2.2e-16", x = 4, y = 1, size = 4, colour = "red") 

####################

### MICHA

head(ProC5)
xxx=cbind(R1_count,C5_count)
xxx
order(xxx, xxx[,2])
order(xxx[,2])
xxx[,order(xxx[,2])]
xxx[order(xxx[,2])]
xxx
xxx[order(xxx[,2]),]
x4=xxx[order(xxx[,2]),]
xs=head(x4,100);cor(xs[,1],xs[,2])
xs
x4=xxx[order(xxx[,2],dec=T),]
x4=xxx[order(xxx[,2],decreasing =T),]
x4
summary(lm(x4[,2]~x4[,1]))
xs=head(x4,100);cor(xs[,1],xs[,2])
xs=head(x4,1000);cor(xs[,1],xs[,2])
xs=head(x4,10);cor(xs[,1],xs[,2])
xs=x4[10:100];cor(xs[,1],xs[,2])
xs=x4[10:100,];cor(xs[,1],xs[,2])
xs=x4[100:1000,];cor(xs[,1],xs[,2])
xs=x4[1000:10000,];cor(xs[,1],xs[,2])
xs=x4[10000:100000,];cor(xs[,1],xs[,2])
xs=x4[10000:nrow(x4),];cor(xs[,1],xs[,2])
plot(xs[,1],xs[,2])
xs[xs[,1]>2000,]
plot(xs[,1],xs[,2],log="xy")
plot(xs[,1],xs[,2])
plot(xs[,1],xs[,2],xlim=c(0,500))
plot(xs[,1],xs[,2],xlim=c(0,250))
plot(xs[,1],xs[,2])
total_counts_median[total_counts_median$ProR1==2714&total_counts_median$ProC5==244,]
xs=x4[100:1000,];cor(xs[,1],xs[,2])
plot(xs)
plot(xs,log="xy")
cor(x4)
cor(x4[,1],x4[,2])
a=1;b=100;cor(x4[a:b,1],x4[a:b,2])
a=100;b=1000;cor(x4[a:b,1],x4[a:b,2])
a=1000;b=10000;cor(x4[a:b,1],x4[a:b,2])
a=10000;b=nrow(x4);cor(x4[a:b,1],x4[a:b,2])
a=10000;b=nrow(x4);cor(x4[a:b,1],x4[a:b,2],log="xy")
a=10000;b=nrow(x4);plot(x4[a:b,1],x4[a:b,2],log="xy")

##################### END MICHA ################################
