#!/usr/bin/env bash
ZIPPED_GTF=$1

LAST_GTF_LINE=`zcat $ZIPPED_GTF | grep 'gene_type "protein_coding"' | awk '$3=="exon"' | wc -l`

zcat $ZIPPED_GTF | grep 'gene_type "protein_coding"' | awk -v FS="\t" -v OFS="\t" '$3=="exon"{

    split($9,source,"\"");
    print $1,$4,$5,source[2],source[4]}' | sort -k4,4 -k5,5 | awk -v LastLine=$LAST_GTF_LINE -v OFS="\t" -v FS="\t" '{

    if(NR==1){

	    gene = $4;
        tx= $5;
        l = $3-$2+1;
        next;
    }

    if(tx==$5){
        l+= $3-$2+1;

    }else{

        print gene,l,tx;
        gene = $4;
        tx = $5;
        l = $3-$2+1;
    }

    if(NR==LastLine){
        print gene,l,tx;
    }
}'
