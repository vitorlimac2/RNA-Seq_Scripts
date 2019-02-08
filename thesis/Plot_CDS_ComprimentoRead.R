### ler arquivo de entrada de ProC/R

setwd("/media/vitor/Seagate Expansion Drive/Thesis/")

### Input file example
## READ_ID  READ_LENGHT CDS_LENGTH  GC_GENE
## SRR2534126.10000274 53 617 43.7908
## SRR2534126.10000276 93 543 45.1673

f <- read.table("ProR1.ReadLength.MedianLenght.GC", header=F)[,2:4]

colnames(f) <- c("length","cds","gc")

head(f)
### criar os breaks iguais

my_breaks <- c(20,50,100,150,200,250)


## plot de 5 vetores para CDS

cds20 <- f[f$length < 50,]$cds
cds50 <- f[f$length >= 50 & f$length < 100,]$cds
cds100 <- f[f$length >= 100 & f$length < 150,]$cds
cds150 <- f[f$length >= 150 & f$length < 200,]$cds
cds200 <- f[f$length >= 200 & f$length < 250,]$cds
cds250 <- f[f$length >= 250,]$cds


## plot de 5 vetores para GC

labels <- c("[20,50[","[50,100[","[100,150[","[150,200[","[200,250[","[250,...[")

boxplot(cds20, 
        cds50, 
        cds100, 
        cds150, 
        cds200, 
        cds250, 
        notch=T,
        log = 'y', 
        outline = F,
        ylab="Comprimento mediano de CDS (pb)",
        xlab="Comprimento de reads mapeados (pb)",
        names=labels)

gc20 <- f[f$length < 50,]$gc
gc50 <- f[f$length >= 50 & f$length < 100,]$gc
gc100 <- f[f$length >= 100 & f$length < 150,]$gc
gc150 <- f[f$length >= 150 & f$length < 200,]$gc
gc200 <- f[f$length >= 200 & f$length < 250,]$gc
gc250 <- f[f$length >= 250,]$gc


## plot de 5 vetores para GC
library("vioplot")
boxplot(gc20, 
        gc50, 
        gc100, 
        gc150, 
        gc200, 
        gc250, 
        notch=T, 
        ylab="Conteúdo GC (%)",
        xlab="Comprimento de reads mapeados (pb)",
        names=labels,
        outline = F)

