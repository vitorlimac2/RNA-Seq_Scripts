samtools view ProC1Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC1.MappedReadInfo

samtools view ProC2Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC2.MappedReadInfo

samtools view ProC3Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC3.MappedReadInfo

samtools view ProC4Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC4.MappedReadInfo

samtools view ProC5Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC5.MappedReadInfo

samtools view ProC6Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC6.MappedReadInfo

samtools view ProC7Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC7.MappedReadInfo

samtools view ProC8Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC8.MappedReadInfo

samtools view ProC9Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProC9.MappedReadInfo

samtools view ProR1Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProR1.MappedReadInfo

samtools view ProR2Aligned.out.bam | awk -v FS="\t" -v OFS="\t" '{read_length[$1]=length($10); if($3!="*"){read_status[$1]=1; split($12,nh,":"); hit_number[$1]=nh[3]}else{if($3!=1){read_status[$1]=0}}} END {for(i in read_length){print i, read_length[i], hit_number[i],read_status[i]}}' > ProR2.MappedReadInfo

