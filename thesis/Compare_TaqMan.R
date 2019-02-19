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

taqman <- read.table(taqman_path, header=T)

raw_file$ProC_2 <- log10(raw_file$ProC_2 + 0.00000001)
raw_file$ProR <- log10(raw_file$ProR  + 0.00000001)
raw_file$HiSeq <- log10(raw_file$HiSeq + 00000001)

norm_file$ProC_2 <- log10(norm_file$ProC_2  + 0.00000001)
norm_file$ProR <- log10(norm_file$ProR  + 0.00000001)
norm_file$HiSeq <- log10(norm_file$HiSeq + 00000001)

taqman$taqman <- log10(taqman$taqman  + 0.00000001)

raw <- merge(raw_file, taqman, by = "gene", sort = T)
norm <- merge(norm_file, taqman, by = "gene", sort = T)

l1 <- lm(raw$ProR ~ norm$taqman); rsq(l1)


boxplot(raw$ProC_2, norm$ProC_2, raw$ProR, norm$ProR, raw$HiSeq, norm$taqman, outline = F)

c1 <- cor.test(raw$ProC_2, raw$taqman, method="spearman")
c2 <- cor.test(norm$ProC_2, norm$taqman, method="spearman")

c3 <- cor.test(raw$ProR, raw$taqman, method="spearman")
c4 <- cor.test(norm$ProR, norm$taqman, method="spearman")

cor.test(log2(raw$HiSeq + 0.00000001), raw$taqman, method="spearman")
