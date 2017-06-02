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
}' | sort -k1,1 -k2,2n > gene_length_tx

LAST_INPUT_LINE=`cat gene_length_tx | wc -l`

cat gene_length_tx | awk -v FS="\t" -v OFS="\t" -v LastLine=$LAST_INPUT_LINE '
function round(x,   ival, aval, fraction)
{
   ival = int(x)    # integer part, int() truncates

   # see if fractional part
   if (ival == x)   # no fraction
      return ival   # ensure no decimals

   if (x < 0) {
      aval = -x     # absolute value
      ival = int(aval)
      fraction = aval - ival
      if (fraction >= .5)
         return int(x) - 1   # -2.5 --> -3
      else
         return int(x)       # -2.3 --> -2
   } else {
      fraction = x - ival
      if (fraction >= .5)
         return ival + 1
      else
         return ival
   }
}
{
    if(NR==1){
       gene = $1;
       l_len[$2"\t"$3] = $2;
       next;
    }

    if(gene==$1){
        l_tx[$2] = $3;
        l_len[$3] = $2;

    }else{

        middle=asort(l_len,dest)/2;

        if(length(l_len)%2==0){
            print gene, dest[middle], l_tx[dest[middle]];
        }else{
            print gene, dest[round(middle)], l_tx[dest[round(middle)]];
        }

        delete l_tx;
        delete l_len;
        delete dest;
        gene = $1;
        l_tx[$2] = $3;
        l_len[$3] = $2;

    }

    if(NR==LastLine){
        middle=asort(l_len,dest)/2;
        if(length(l_len)%2==0){
            print gene, dest[middle], l_tx[dest[middle]];
        }else{
            print gene, dest[round(middle)], l_tx[dest[round(middle)]];

        }

    }
}'