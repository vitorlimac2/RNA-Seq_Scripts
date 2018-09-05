### Plot Mapping Rate By Read Length

ProR_2 <- c(85.4, 92.1, 88.8, 80.8, 72.2, 63.9)

ProC_2 <- c(83.3, 82.2, 84.6, 89.9, 86.4, 80.3)

ProC_8 <- c(76.7, 72 , 73.8, 78.4, 71.9, 74.5)

ranges <- c("[20,50[", "[50,100[", "[100,150[",
            "[150,200[", "[200,250[","[250,...[")

mappings <- data.frame(ProR_2, ProC_2, ProC_8, ranges)


mappings$ranges <- factor(mappings$ranges, levels=ranges)

library("ggplot2")
library("tidyr")

p <- mappings %>% gather(Replicate, Percentage, ProR_2, ProC_2, ProC_8) %>% ggplot(aes(x=ranges, y=Percentage, colour=Replicate, group=Replicate)) + geom_point() + geom_line() + theme_bw()

p <- p + labs(x="Comprimento de read (pb)", y="%")

p + geom_text(aes(label=Percentage, hjust=1, vjust=2))
