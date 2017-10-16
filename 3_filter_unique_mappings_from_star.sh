## BAM file of STAR
a=$1
xpath=${a%/*}
xbase=${a##*/}
xfext=${xbase##*.}
xpref=${xbase%.*}

samtools view $a | awk -v FS="\t" -v OFS="\t" '$5==255' | sort -k1,1 > $xpref.unique.usort.sam
