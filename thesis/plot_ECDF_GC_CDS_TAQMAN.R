setwd("/media/vitor/Seagate Expansion Drive/Thesis/")


low_raw_proc <- read.table("ProC_rna2.LowFragmentation.CDS_median.GC")
high_raw_proc <- read.table("ProC_rna2.HighFragmentation.CDS_median.GC")
low_raw_pror <- read.table("ProR.LowFragmentation.CDS_median.GC")
high_raw_pror <- read.table("ProR.HighFragmentation.CDS_median.GC")
background_raw <- read.table("Homo_sapiens.GRCh38.CDS_median.GC")


low_pror <- low_raw_pror[sample(nrow(low_raw_pror), .25*nrow(low_raw_pror)),]$V4
high_pror <- high_raw_pror[sample(nrow(high_raw_pror), .25*nrow(high_raw_pror)),]$V4
low_proc <- low_raw_proc[sample(nrow(low_raw_proc), .25*nrow(low_raw_proc)),]$V4
high_proc <- high_raw_proc[sample(nrow(high_raw_proc), .25*nrow(high_raw_proc)),]$V4
background <- background_raw[sample(nrow(background_raw), .75*nrow(background_raw)),]$V4

plot(main="Genes com perfil de alta fragmentação", 
     ecdf(high_proc), 
     col = 'blue', 
     lwd=7, 
     ylab = "Genes (%)", 
     xlab="Conteúdo GC (%)",
     yaxt='n')

axis(2, at=c(0, .2, .4, .6, .8, 1), labels = c(0, 20, 40, 60, 80, 100))

plot(ecdf(high_pror), add=TRUE, col = 'darkblue', lwd=7)
plot(ecdf(background), add=TRUE, col = 'red', lwd=7)

legend(
       "topleft", 
       legend = c("Background",
                             "ProR",
                             "ProC"),
       col=c("red","darkblue","blue"), 
       lty=c(1,1,1)
       , bty = 'n', lwd=7)


plot(main="Genes com perfil de baixa fragmentação",
     ecdf(low_proc), 
     col='blue', 
     lwd=1, 
     ylab = "Genes (%)", 
     xlab="Conteúdo GC (%)",
     yaxt='n')

## 0, .2, .4, .6, .8, 1
axis(2, at=c(0, .2, .4, .6, .8, 1), labels = c(0, 20, 40, 60, 80, 100))

plot(ecdf(low_pror), add=TRUE, col = 'darkblue', lwd=1)
plot(ecdf(background), add=TRUE, col = 'red', lwd=7)

legend(
       "topleft", 
       legend = c("Background",
                  "ProR",
                  "ProC"),
       col=c("red","darkblue","blue"), 
       lty=c(1,1,1)
       , bty = 'n', lwd=7)

################################################################################

low_pror <- low_raw_pror[sample(nrow(low_raw_pror), .25*nrow(low_raw_pror)),]$V3
high_pror <- high_raw_pror[sample(nrow(high_raw_pror), .25*nrow(high_raw_pror)),]$V3
low_proc <- low_raw_proc[sample(nrow(low_raw_proc), .25*nrow(low_raw_proc)),]$V3
high_proc <- high_raw_proc[sample(nrow(high_raw_proc), .25*nrow(high_raw_proc)),]$V3
background <- background_raw[sample(nrow(background_raw), .75*nrow(background_raw)),]$V3

plot(main="Genes com perfil de alta fragmentação", 
     ecdf(log10(high_proc)), 
     col = 'blue', 
     lwd=7, 
     ylab = "Genes (%)", 
     xlab="CDS log10(pb)",
     yaxt='n')

axis(2, at=c(0, .2, .4, .6, .8, 1), labels = c(0, 20, 40, 60, 80, 100))

plot(ecdf(log10(high_pror)), add=TRUE, col = 'darkblue', lwd=7)
plot(ecdf(log10(background)), add=TRUE, col = 'red', lwd=7)

legend(
  "topleft", 
  legend = c("Background",
             "ProR",
             "ProC"),
  col=c("red","darkblue","blue"), 
  lty=c(1,1,1)
  , bty = 'n', lwd=7)


plot(main="Genes com perfil de baixa fragmentação",
     ecdf(log10(low_proc)), 
     col='blue', 
     lwd=1, 
     ylab = "Genes (%)", 
     xlab="CDS log10(pb)",
     yaxt='n')

## 0, .2, .4, .6, .8, 1
axis(2, at=c(0, .2, .4, .6, .8, 1), labels = c(0, 20, 40, 60, 80, 100))

plot(ecdf(log10(low_pror)), add=TRUE,col = 'darkblue', lwd=1)
plot(ecdf(log10(background)), add=TRUE, col = 'red', lwd=7)

legend(
  "topleft", 
  legend = c("Background",
             "ProR",
             "ProC"),
  col=c("red","darkblue","blue"), 
  lty=c(1,1,1)
  , bty = 'n', lwd=7)



######################################3 MAQC qRT-PCR genes

maqc <- read.table("MAQC_ENSEMBL_CDS")
background_raw <- read.table("Homo_sapiens.GRCh38.89_CDS_median_length")

head(maqc)
head(background_raw)

maqc_cds <- log10(maqc$V3)
background_cds <- log10(background_raw$V2)


my_breaks <- pretty(range(c(maqc_cds, background_cds)), n=8)

h_1 <- hist(maqc_cds, breaks=my_breaks, plot=FALSE)
h_2 <- hist(background_cds, breaks=my_breaks, plot=FALSE)


max_y <- max(h_1$counts*100/length(maqc_cds), 
             h_2$counts*100/length(background_cds))

plot(h_1$mids, h_1$counts*100/length(maqc_cds), 
     type="l", 
     ylab = "%", 
     xlab = "Comprimento de CDS  (nt)", 
     col="blue", 
     ylim=c(0,max_y),
     lwd=6) 

  lines(h_2$mids, h_2$counts*100/length(background_cds), col="red")



legend("topleft", legend = c("Background",
                             "Baixa Fragmentação",
                             "Alta Fragmentação"),
       col=c("gray22","black","blue"), lty=c(4,1,1), bty = 'n', lwd=7)


###########################################################################################

maqc <- read.table("MAQC_ENSEMBL_CDS_GC")
background_raw <- read.table("Homo_sapiens.GRCh38.89.exon.Tx.ProteinCodingGenesAndTx.GC")

head(maqc)
head(background_raw)

maqc_gc <- log10(maqc$V4)
background_gc <- log10(background_raw$V2)


my_breaks <- pretty(range(c(maqc_gc, background_gc)), n=5)

h_1 <- hist(maqc_gc, breaks=my_breaks, plot=FALSE)
h_2 <- hist(background_gc, breaks=my_breaks, plot=FALSE)


max_y <- max(h_1$counts*100/length(maqc_gc), 
             h_2$counts*100/length(background_gc))

plot(h_1$mids, h_1$counts*100/length(maqc_gc), 
     type="l", 
     ylab = "%", 
     xlab = "Comprimento de CDS  (nt)", 
     col="blue", 
     ylim=c(0,max_y),
     lwd=6) 

lines(h_2$mids, h_2$counts*100/length(background_gc), col="red")



legend("topleft", legend = c("Background",
                             "Baixa Fragmentação",
                             "Alta Fragmentação"),
       col=c("gray22","black","blue"), lty=c(4,1,1), bty = 'n', lwd=7)

#############################################################################################################################
#############################################################################################################################
### PLOT DEG MUS MUSCULUS

setwd("/media/vitor/Seagate Expansion Drive/Thesis/DEG_ILB/")


GC <- 4
CDS <- 3

### DEG ON BOTH QUANTIFICATIONS

both  <- read.table("Genes_UP-DEG_Ambos.CdsMedianLength.GC")

#### OLD DEG ON NOT-ADJUSTED QUANT
deg_nadj <- read.table("Genes_UP-DEG_apenas_no_Raw.CdsMedianLength.GC")


#### NEW DEG ON ADJUSTED QUANT.
deg_adj <- read.table("Genes_UP-DEG_apenas_no_Norm.CdsMedianLength.GC")

f <- read.table("Mus_musculus.NCBIM37.67.protein_coding.median_CDS.gc_CDS")

#background <- both[,CDS]

par(mar= c(5,4,2,2))

boxplot(deg_adj[,CDS], 
        deg_nadj[,CDS], 
        background, 
        outline = F, 
        log='y' ,
        notch = F,
        ylab="Comprimento de CDS (pb)",
        names=c("Novos DEG",
                 "DEG somente \nna quant. não-ajustada",
                 "Background"))



#background <- both[sample(nrow(both), 100),CDS]


background <- f[sample(nrow(f), 500),]

par(mar= c(5,4,2,2))

boxplot(deg_adj[,GC], 
        deg_nadj[,GC], 
        background, 
        outline = F, 
        notch = F,
        ylab="Conteúdo GC (%)",
        names=c("Novos DEG",
                "DEG somente \nna quant. não-ajustada",
                "Background"))
#######

my_breaks <- pretty(range(c(deg_adj[,GC],deg_nadj[,GC], both[,GC]), background[,GC]), n = 4)

h1 <- hist(deg_adj[,GC], breaks=my_breaks, plot = F)
h2 <- hist(deg_nadj[,GC], breaks=my_breaks, plot = F)
h3 <- hist(both[,GC], breaks=my_breaks, plot = F)
h4 <- hist(background[,GC], breaks=my_breaks, plot = F)

h1$mids;
h1$counts*100/length(deg_adj[,GC])
h2$counts*100/length(deg_nadj[,GC])
h3$counts*100/length(both[,GC])
h4$counts*100/length(background[,GC])

my_breaks <- pretty(range(c(deg_adj[,CDS],deg_nadj[,CDS], both[,CDS], background[,CDS])), n =5)

h1 <- hist(deg_adj[,CDS], breaks=my_breaks, plot = F)
h2 <- hist(deg_nadj[,CDS], breaks=my_breaks, plot = F)
h3 <- hist(both[,CDS], breaks=my_breaks, plot = F)
h4 <- hist(background[,CDS], breaks=my_breaks, plot = F)

plot(ecdf)


h1$mids;
h1$counts*100/length(deg_adj[,CDS])
h2$counts*100/length(deg_nadj[,CDS])
h3$counts*100/length(both[,CDS])
h4$counts*100/length(background[,CDS])

#########################################################################################################################
############# GENES COM MAIOR MUDANÇA DE FOLD-CHANGE

setwd("/media/vitor/Seagate Expansion Drive/Thesis/DEG_ILB/")


GC <- 4
CDS <- 3

### DEG ON BOTH QUANTIFICATIONS

f <- read.table("Mus_musculus.NCBIM37.67.protein_coding.median_CDS.gc_CDS")

both  <- read.table("Genes_UP-DEG_Ambos.CdsMedianLength.GC")

deg_altoFC <- read.table("Genes_Fold-Change_maior_que_1.CdsMedianLength.GC")

GC <- 4
CDS <- 3

#my_breaks <- pretty(range(c(deg_altoFC[,GC],f[,GC]), n = 5))
my_breaks <- pretty(range(c(deg_altoFC[,GC],both[,GC]), n = 5))


h1 <- hist(deg_altoFC[,GC], breaks=my_breaks, plot = F)
h2 <- hist(f[,GC], breaks=my_breaks, plot = F)
h3 <- hist(both[,GC], breaks=my_breaks, plot = F)


h1$mids;
h1$counts*100/length(deg_altoFC[,GC])

h2$counts*100/length(f[,GC])

h3$counts*100/length(both[,GC])

#####################################################

my_breaks <- pretty(range(c(deg_altoFC[,CDS],both[,CDS]), n = 30))


h1 <- hist(deg_altoFC[,CDS], breaks=my_breaks, plot = F)
h3 <- hist(both[,CDS], breaks=my_breaks, plot = F)

boxplot(deg_adj[,CDS], deg_nadj[,CDS],both[,CDS], outline = F, ylab="Comprimento CDS (pb)", names = c("Novos GDEs","Somente quant.\nnão-ajustada",  "Ambas quant."), 
        las=1,cex.axis=0.8)

h1$mids;
h1$counts*100/length(deg_altoFC[,CDS])
h3$counts*100/length(both[,CDS])
