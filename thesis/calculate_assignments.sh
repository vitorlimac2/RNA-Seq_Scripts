
REP=$1

join "$REP".MappedReadInfo "$REP"Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$3==1 && $5!="NA"{print $1, $2}' | sort | uniq  > temp

cat temp | awk '$2 < 50' | wc -l
cat temp | awk '$2 >= 50 && $2 < 100' | wc -l
cat temp | awk '$2 >= 100 && $2 < 150' | wc -l
cat temp | awk '$2 >= 150 && $2 < 200' | wc -l
cat temp | awk '$2 >= 200 && $2 < 250' | wc -l
cat temp | awk '$2 >= 250' | wc -l

