sort -k1,1 ProC1.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC1.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC2.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC2.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC3.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC3.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC4.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC4.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC5.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC5.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC6.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC6.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC7.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC7.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC8.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC8.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProC9.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProC9.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProR1.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProR1.MappedReadInfo.ReadCountNormalizedByRL

sort -k1,1 ProR2.MappedReadInfo.ReadCountNormalizedByRL > temp
mv temp ProR2.MappedReadInfo.ReadCountNormalizedByRL

join ProC1.MappedReadInfo.ReadCountNormalizedByRL ProC1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC1.normalized.Frac

join ProC1.MappedReadInfo.ReadCountNormalizedByRL ProC1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC1.normalized.NoFrac

join ProC1.MappedReadInfo.ReadCountNormalizedByRL ProC1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC1.NoNormalized.NoFrac

join ProC1.MappedReadInfo.ReadCountNormalizedByRL ProC1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC1.NoNormalized.Frac

join ProC2.MappedReadInfo.ReadCountNormalizedByRL ProC2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC2.normalized.Frac

join ProC2.MappedReadInfo.ReadCountNormalizedByRL ProC2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC2.normalized.NoFrac

join ProC2.MappedReadInfo.ReadCountNormalizedByRL ProC2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC2.NoNormalized.NoFrac

join ProC2.MappedReadInfo.ReadCountNormalizedByRL ProC2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC2.NoNormalized.Frac

join ProC3.MappedReadInfo.ReadCountNormalizedByRL ProC3Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC3.normalized.Frac

join ProC3.MappedReadInfo.ReadCountNormalizedByRL ProC3Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC3.normalized.NoFrac

join ProC3.MappedReadInfo.ReadCountNormalizedByRL ProC3Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC3.NoNormalized.NoFrac

join ProC3.MappedReadInfo.ReadCountNormalizedByRL ProC3Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC3.NoNormalized.Frac

join ProC4.MappedReadInfo.ReadCountNormalizedByRL ProC4Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC4.normalized.Frac

join ProC4.MappedReadInfo.ReadCountNormalizedByRL ProC4Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC4.normalized.NoFrac

join ProC4.MappedReadInfo.ReadCountNormalizedByRL ProC4Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC4.NoNormalized.NoFrac

join ProC4.MappedReadInfo.ReadCountNormalizedByRL ProC4Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC4.NoNormalized.Frac

join ProC5.MappedReadInfo.ReadCountNormalizedByRL ProC5Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC5.normalized.Frac

join ProC5.MappedReadInfo.ReadCountNormalizedByRL ProC5Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC5.normalized.NoFrac

join ProC5.MappedReadInfo.ReadCountNormalizedByRL ProC5Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC5.NoNormalized.NoFrac

join ProC5.MappedReadInfo.ReadCountNormalizedByRL ProC5Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC5.NoNormalized.Frac

join ProC6.MappedReadInfo.ReadCountNormalizedByRL ProC6Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC6.normalized.Frac

join ProC6.MappedReadInfo.ReadCountNormalizedByRL ProC6Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC6.normalized.NoFrac

join ProC6.MappedReadInfo.ReadCountNormalizedByRL ProC6Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC6.NoNormalized.NoFrac

join ProC6.MappedReadInfo.ReadCountNormalizedByRL ProC6Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC6.NoNormalized.Frac

join ProC7.MappedReadInfo.ReadCountNormalizedByRL ProC7Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC7.normalized.Frac

join ProC7.MappedReadInfo.ReadCountNormalizedByRL ProC7Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC7.normalized.NoFrac

join ProC7.MappedReadInfo.ReadCountNormalizedByRL ProC7Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC7.NoNormalized.NoFrac

join ProC7.MappedReadInfo.ReadCountNormalizedByRL ProC7Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC7.NoNormalized.Frac

join ProC8.MappedReadInfo.ReadCountNormalizedByRL ProC8Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC8.normalized.Frac

join ProC8.MappedReadInfo.ReadCountNormalizedByRL ProC8Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC8.normalized.NoFrac

join ProC8.MappedReadInfo.ReadCountNormalizedByRL ProC8Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC8.NoNormalized.NoFrac

join ProC8.MappedReadInfo.ReadCountNormalizedByRL ProC8Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC8.NoNormalized.Frac

join ProC9.MappedReadInfo.ReadCountNormalizedByRL ProC9Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC9.normalized.Frac

join ProC9.MappedReadInfo.ReadCountNormalizedByRL ProC9Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC9.normalized.NoFrac

join ProC9.MappedReadInfo.ReadCountNormalizedByRL ProC9Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC9.NoNormalized.NoFrac

join ProC9.MappedReadInfo.ReadCountNormalizedByRL ProC9Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProC9.NoNormalized.Frac

join ProR1.MappedReadInfo.ReadCountNormalizedByRL ProR1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR1.normalized.Frac

join ProR1.MappedReadInfo.ReadCountNormalizedByRL ProR1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR1.normalized.NoFrac

join ProR1.MappedReadInfo.ReadCountNormalizedByRL ProR1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR1.NoNormalized.NoFrac

join ProR1.MappedReadInfo.ReadCountNormalizedByRL ProR1Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR1.NoNormalized.Frac

join ProR2.MappedReadInfo.ReadCountNormalizedByRL ProR2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR2.normalized.Frac

join ProR2.MappedReadInfo.ReadCountNormalizedByRL ProR2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=$5} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR2.normalized.NoFrac

join ProR2.MappedReadInfo.ReadCountNormalizedByRL ProR2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR2.NoNormalized.NoFrac

join ProR2.MappedReadInfo.ReadCountNormalizedByRL ProR2Aligned.out.bam.featureCounts.Unique.Sorted.ReadAssignmentToGene | awk '$7!=0{gene_count[$6]+=1/$7} END {for(i in gene_count){print i, gene_count[i]}}' | sort -k1,1 > ProR2.NoNormalized.Frac

