setwd("/media/vitor/Seagate Expansion Drive/Thesis/")

norm_path = "NormalizedAllMean.tsv" 
raw_path = "RawAllMean.tsv"
taqman_path <- "ENSEMBL_MEAN_EXPRESSION_ID_SYMBOL"

## read raw file
raw_file <- read.table(raw_path, 
                       header = T)
## read norm file
norm_file <- read.table(norm_path,
                        header=T)

raw_file$ProC_2 <- raw_file$ProC_2*1000000/sum(raw_file$ProC_2)
norm_file$ProC_2 <- norm_file$ProC_2*1000000/sum(norm_file$ProC_2)

raw_file$ProR <- raw_file$ProR*1000000/sum(raw_file$ProR)
norm_file$ProR <- norm_file$ProR*1000000/sum(norm_file$ProR)

raw_file$ProC_2 <- log10(raw_file$ProC_2 + 0.00000001)
norm_file$ProC_2 <- log10(norm_file$ProC_2  + 0.00000001)

raw_file$ProR <- log10(raw_file$ProR  + 0.00000001)
norm_file$ProR <- log10(norm_file$ProR  + 0.00000001)

taqman <- read.table(taqman_path, header=T)


x <- merge(raw_file, taqman, by = "gene", sort = T)

z <- merge(norm_file, taqman, by = "gene", sort = T)

c1 <- cor.test(x$ProC_2, x$taqman, method="spearman")
c2 <- cor.test(z$ProC_2, x$taqman, method="spearman")

c3 <- cor.test(x$ProR, x$taqman, method="spearman")
c4 <- cor.test(z$ProR, x$taqman, method="spearman")

plot(x$ProR, x$taqman)

head(x)


c1
c2
c3
c4

l1 <- lm(x$ProC_2 ~ x$taqman)
summary(l1)$r.squared

l2 <- lm(z$ProC_2 ~ x$taqman)
summary(l2)$r.squared
