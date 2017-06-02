# Receive a input with three columns

# Count_of_occurrence   Read_name	flag

## Flag 4 = unmapped
## Flag 16 = mapped to reverse strand
## Flag 0 = mapped to forward strand


INPUT_SAM=$1

ORDERED_READS = `awk '{print $1,$2}' $INPUT_SAM | grep -v "^@SQ" | sort | uniq -c`

awk -v unique=0 -v multi=0 -v unmapped=0 -v mapped=0 '{

if($3 == 4){
	unmapped=unmapped+1

}else{

read_names[$2] = read_names[$2] + 1 

}

}

END

{

for(i in read_names){
	mapped=mapped+1;
	if(read_names[i]==1){
		unique=unique+1
	}else{
		multi=multi+1
	}

}

print "Total Reads = "mapped+unmapped
print "\tUnique mapped reads = "unique" "unique*100/(mapped+unmapped)
print "\tMulti-mapped reads = "multi" "multi*100/(mapped+unmapped)
print "\tUnmapped reads = "unmapped" "unmapped*100/(mapped+unmapped)

}' $ORDERED_READS
