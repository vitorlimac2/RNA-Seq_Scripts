
total_counts <- read.table("ProR1_ProC5.sorted.total")

plot(log10(total_counts$V2+1),log10(total_counts$V3+1), ylab = "ProC5", xlab = "ProR1", main = "Counts (log10)", pch = 16, cex = 0.5, col = "blue")

model1 <- lm (log10(total_counts$V3+1)  ~ log10(total_counts$V2+1))

abline(model1)

summary(model1)

text(4,1,"R² = 0.9; p-value < 2.2e-16")


plot(log10(total_counts$V2+1),log10(total_counts$V3+1), ylab = "ProC5", xlab = "ProR1", main = "Counts (log10)", pch = 16, cex = 0.5, col = "black")
abline(model1, lwd=2)

res <- signif(residuals(model1), 5)
pre <- predict(model1)

segments(log10(total_counts$V2+1), log10(total_counts$V3+1), log10(total_counts$V2+1), pre, col="red")

library(calibrate)
textxy(log10(total_counts$V2+1), log10(total_counts$V3+1), rownames(total_counts), cx=0.7)


#####################################################################

total_counts_median <- read.table("ProR1_ProC5.sorted.total.median_length")
model1 <- lm (log10(total_counts_median$V3+1)  ~ log10(total_counts_median$V2+1))
library(ggplot2)

colnames(total_counts_median) <- c("gene", "ProR1", "ProC5", "Length")

attach(total_counts_median)

p <- ggplot(total_counts_median, aes(x = log10(ProR1+1), 
                                y = log10(ProC5+1), 
                                color = "blue",
                                size = log10(Length))) +
  geom_point()

p + geom_abline(intercept = 0.09399136, slope = 1.02267463) + annotate("text", label = "R² = 0.9; p-value < 2.2e-16", x = 4, y = 1, size = 4, colour = "red")

ggplot(total_counts_median, aes(x = log10(ProR1+1), 
                                y = log10(ProC5+1), 
                                size = log10(Length))) +
  geom_point()

plot(log10(total_counts$V2+1),log10(total_counts$V3+1), ylab = "ProC5", xlab = "ProR1", main = "Counts (log10)", pch = 16, cex = 0.5, col = "black")
abline(model1)
res <- signif(residuals(model1), 5)
pre <- predict(model1)
library(calibrate)
textxy(log10(total_counts$V2+1), log10(total_counts$V3+1), total_counts_median$V4, cx=0.7)

xx=10;cor(head(sort(log10(ProC5+1),dec=T),xx), head(sort(log10(ProR1+1),dec=T),xx))
head(which(is.na(total_counts_median$V3)))

head(total_counts_median$V3)

####################

### MICHA

head(ProC5)
xxx=cbind(ProR1,ProC5)
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
