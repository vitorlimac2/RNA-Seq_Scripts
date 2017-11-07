#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=1){
  message("USAGE:\nRscript --vanilla KS_test_Dist_Frag_Qvalue.R ks_pvalues_output_file.\nOPTIONS:\n\tks_pvalues_output_file: Output file of KS_test_Dist_Frag.R with only p-values on columns", call.=FALSE)
  stop("Missing options.")
}


inputFile <- args[1]
pvalues_table <- read.table(inputFile, header=F, row.names=1)
colnames(pvalues_table) <- c("less","two_sided","greater")
nr <- nrow(pvalues_table)
qvalues_table <- pvalues_table
qvalues_table$less <- p.adjust(pvalues_table$less, method="fdr", n=nr)
qvalues_table$two_sided <- p.adjust(pvalues_table$two_sided, method="fdr", n=nr)
qvalues_table$greater <- p.adjust(pvalues_table$greater, method="fdr", n=nr)
output_file <- paste(inputFile,".Qvalues",sep="")
write.table(qvalues_table,output_file,sep="\t")

