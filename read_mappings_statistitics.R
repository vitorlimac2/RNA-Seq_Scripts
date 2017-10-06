setwd("~/PRJ.SRP064142/QUANTIFY_REPLICATES")
file_suffix <- "_read_mps_trim_raw"

my_replicates <- c("ProC1",
                                 "ProC2",
                                 "ProC3",
                                 "ProC4",
                                 "ProC5",
                                 "ProC6",
                                 "ProC7",
                                 "ProC8",
                                 "ProC9",
                                 "ProR1",
                                 "ProR2")

for(i in my_replicates){
  
  rep_preffix <- i
  
  
  my_file_name <- paste(rep_preffix,file_suffix,sep="")
  
  my_rep <- read.table(my_file_name, 
                       row.names = 1, 
                       colClasses = c("character","numeric", "numeric","numeric"))
  
  ## hist of total mapped; unique; multi
  
  mapped_reads <- my_rep[my_rep$V2 > 0,]$V3
  unique_reads <- my_rep[my_rep$V2 == 1,]$V3
  mm_reads <- my_rep[my_rep$V2 > 1,]$V3
  un_reads <- my_rep[my_rep$V2 == 0,]$V3
  
  my_breaks <- pretty(range(c(my_rep$V4,mapped_reads, unique_reads, mm_reads, un_reads)), n=20)
  
  total_reads.h <-hist(my_rep$V3, breaks=my_breaks, plot=FALSE)
  mapped_reads.h <- hist(mapped_reads, breaks=my_breaks, plot=FALSE)
  unique_reads.h <- hist(unique_reads, breaks=my_breaks, plot=FALSE)
  mm_reads.h <- hist(mm_reads, breaks=my_breaks, plot=FALSE)
  un_reads.h <- hist(un_reads, breaks=my_breaks, plot=FALSE)
  
  my_classes <- c("All", "Mapped", "Unique", "Multi-mapped", "Unmapped")
  
  #max_y <- max(total_reads.h$counts*100/nrow(my_rep),
   #            mapped_reads.h$counts*100/nrow(my_rep),
    #           unique_reads.h$counts*100/nrow(my_rep),
     #          mm_reads.h$counts*100/nrow(my_rep),
      #         un_reads.h$counts*100/nrow(my_rep))
  
  max_y <- 40
  
  #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  
  png(paste(my_file_name,".png",sep=""))
  
  plot(main=rep_preffix,total_reads.h$mids, total_reads.h$counts*100/nrow(my_rep), 
       type="l", 
       ylim=c(0,max_y),
       xlab="Read Length (bp)",
       ylab="%")
  
  lines(mapped_reads.h$mids, mapped_reads.h$counts*100/nrow(my_rep), col="blue")
  lines(unique_reads.h$mids, unique_reads.h$counts*100/nrow(my_rep), col="darkgreen")
  lines(mm_reads.h$mids, mm_reads.h$counts*100/nrow(my_rep), col="purple")
  lines(un_reads.h$mids, un_reads.h$counts*100/nrow(my_rep), col="orange")
  
  legend(280,15, cex=0.8, my_classes, fill=c("black",
                                             "blue",
                                             "darkgreen",
                                             "purple",
                                             "orange",
                                             "yellow"), horiz=FALSE)
  
  dev.off()
  
  rm("my_rep", "mapped_reads", "un_reads", "mm_reads","unique_reads")
  
}