#!/bin/bash

FILE=$1

#awk '{ if(NR%3==0) print $0 }' $FILE  | awk '{ print $4 }' \
#    | awk '{ for(i=0;i<NF;i++) {num=$0 ; sum += num; sumsqr += num^2}} END \
#           {printf "Avg score: %f , std: %f out of %d \n", sum/NR, sqrt((sumsqr/NR)-((sum/NR)^2)), NR}'
awk '{ alignment=$0; n=split(alignment,array,"\t"); if (n > 2) {query = array[1]; score = array[3]; if(score != 0) {printf "%d %f \n", length(query), score }}}' $FILE \
    | awk ' {split($0, array, " "); qlen = array[1]; num=array[2]/qlen; counter += 1; sum += num; sumsqr += num^2; sum_bps+=qlen} END \
           {printf "Avg score: %f , std: %f out of %d reads of total %d bps\n", sum/counter, sqrt((sumsqr/counter)-((sum/counter)^2)), counter, sum_bps}'
