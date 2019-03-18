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

##############################################################################
### correlation by length

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

MAQC_lengths <- read.table("MAQC_ENSEMBL_CDS")[,2:3]
colnames(MAQC_lengths) <- c("gene", "length")

taqman_with_lengths <- merge(taqman, MAQC_lengths)

raw_vs_taqman <- merge(raw_file, taqman_with_lengths, by = "gene", sort = T)

cor86_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length < 300,]$ProC_2), 
                  log2(raw_vs_taqman[raw_vs_taqman$length < 300,]$taqman), method="spearman") 
cor300_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 300 & raw_vs_taqman$length < 500,]$ProC_2), 
                   log2(raw_vs_taqman[raw_vs_taqman$length >= 300 & raw_vs_taqman$length < 500,]$taqman), method="spearman")
cor500_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 500 & raw_vs_taqman$length < 750,]$ProC_2), 
                   log2(raw_vs_taqman[raw_vs_taqman$length >= 500 & raw_vs_taqman$length < 750,]$taqman), method="spearman")
cor750_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >=750 & raw_vs_taqman$length < 1000,]$ProC_2), 
                   log2(raw_vs_taqman[raw_vs_taqman$length >=750 & raw_vs_taqman$length < 1000,]$taqman), method="spearman")
cor1000_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 1000 & raw_vs_taqman$length < 1500,]$ProC_2), 
                    log2(raw_vs_taqman[raw_vs_taqman$length >= 1000 & raw_vs_taqman$length < 1500,]$taqman), method="spearman")
cor1500_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 1500,]$ProC_2), 
                    log2(raw_vs_taqman[raw_vs_taqman$length >= 1500,]$taqman), method="spearman")

norm_vs_taqman <- merge(norm_file, taqman_with_lengths, by="gene", sort = T)

cor86_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length < 300,]$ProC_2), 
                    log2(norm_vs_taqman[norm_vs_taqman$length < 300,]$taqman), method="spearman") 
cor300_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 300 & norm_vs_taqman$length < 500,]$ProC_2), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >= 300 & norm_vs_taqman$length < 500,]$taqman), method="spearman")
cor500_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 500 & norm_vs_taqman$length < 750,]$ProC_2), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >= 500 & norm_vs_taqman$length < 750,]$taqman), method="spearman")
cor750_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >=750 & norm_vs_taqman$length < 1000,]$ProC_2), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >=750 & norm_vs_taqman$length < 1000,]$taqman), method="spearman")
cor1000_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 1000 & norm_vs_taqman$length < 1500,]$ProC_2), 
                      log2(norm_vs_taqman[norm_vs_taqman$length >= 1000 & norm_vs_taqman$length < 1500,]$taqman), method="spearman")
cor1500_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 1500,]$ProC_2), 
                      log2(norm_vs_taqman[norm_vs_taqman$length >= 1500,]$taqman), method="spearman")

cor86_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length < 300,]$ProR), 
                    log2(raw_vs_taqman[raw_vs_taqman$length < 300,]$taqman), method="spearman") 
cor300_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 300 & raw_vs_taqman$length < 500,]$ProR), 
                     log2(raw_vs_taqman[raw_vs_taqman$length >= 300 & raw_vs_taqman$length < 500,]$taqman), method="spearman")
cor500_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 500 & raw_vs_taqman$length < 750,]$ProR), 
                     log2(raw_vs_taqman[raw_vs_taqman$length >= 500 & raw_vs_taqman$length < 750,]$taqman), method="spearman")
cor750_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >=750 & raw_vs_taqman$length < 1000,]$ProR), 
                     log2(raw_vs_taqman[raw_vs_taqman$length >=750 & raw_vs_taqman$length < 1000,]$taqman), method="spearman")
cor1000_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 1000 & raw_vs_taqman$length < 1500,]$ProR), 
                      log2(raw_vs_taqman[raw_vs_taqman$length >= 1000 & raw_vs_taqman$length < 1500,]$taqman), method="spearman")
cor1500_r <- cor.test(log2(raw_vs_taqman[raw_vs_taqman$length >= 1500,]$ProR), 
                      log2(raw_vs_taqman[raw_vs_taqman$length >= 1500,]$taqman), method="spearman")

cor86_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length < 300,]$ProR), 
                    log2(norm_vs_taqman[norm_vs_taqman$length < 300,]$taqman), method="spearman") 
cor300_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 300 & norm_vs_taqman$length < 500,]$ProR), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >= 300 & norm_vs_taqman$length < 500,]$taqman), method="spearman")
cor500_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 500 & norm_vs_taqman$length < 750,]$ProR), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >= 500 & norm_vs_taqman$length < 750,]$taqman), method="spearman")
cor750_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >=750 & norm_vs_taqman$length < 1000,]$ProR), 
                     log2(norm_vs_taqman[norm_vs_taqman$length >=750 & norm_vs_taqman$length < 1000,]$taqman), method="spearman")
cor1000_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 1000 & norm_vs_taqman$length < 1500,]$ProR), 
                      log2(norm_vs_taqman[norm_vs_taqman$length >= 1000 & norm_vs_taqman$length < 1500,]$taqman), method="spearman")
cor1500_n <- cor.test(log2(norm_vs_taqman[norm_vs_taqman$length >= 1500,]$ProR), 
                      log2(norm_vs_taqman[norm_vs_taqman$length >= 1500,]$taqman), method="spearman")
