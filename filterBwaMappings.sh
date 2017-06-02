#!/usr/bin/env bash

INPUT_SAM="/home/vitor/test_sam"

SAM_LINES=`cat $INPUT_SAM | awk '!/^@SQ/&&$2!=4' | wc -l`


cat $INPUT_SAM | awk -v myLastLine=$SAM_LINES '{

    if($2==4){
            print "Remove unmapped reads from input file and sort by read name." > "/dev/stderr"
            exit 1
    }

    if(NR==1){
        read_id=$1
        sam_entry=$0;
        split($12,AS,":");
        alignment_score=AS[3];
        next;
    }

    if($1!=read_id){
        print sam_entry;
        read_id=$1
        sam_entry=$0;
        split($12,AS,":");
        alignment_score=AS[3];
    }else{
        split($12,AS,":");
        if(AS[3]>alignment_score){
            sam_entry=$0
            alignment_score=AS[3]
        }
    }

    if(NR==myLastLine){
        print sam_entry;
    }
}'



