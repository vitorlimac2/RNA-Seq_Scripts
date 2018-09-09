#### Make a boxplot from a single vector 

library("dplyr")
library("ggplot2")

male = data.frame(value=c(127,44,28,83,0,6,78,6,5,213,73,20,214,28,11)) # data from page 66

ggplot(data = male, aes(x = "", y = male$value)) + geom_boxplot() + coord_cartesian(ylim = c(0, 150)) # I set the y axis scale so the plot looks better.

f <- read.table("/home/vitor/Downloads/Tcruzi_CLBrener_Paralogy_TriTrypDB28_groups_size")$V1

f <- data.frame(value=f[f>3])

ggplot(data = f, aes(x = "", y = log10(f$value))) + geom_boxplot() + coord_cartesian() # I set the y axis scale so the plot looks better.

boxplot(log10(f), outline = F)
f

f
