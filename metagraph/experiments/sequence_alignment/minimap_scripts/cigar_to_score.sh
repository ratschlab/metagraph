#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./cigar_to_score <sam filename>" # <match_score> <mismatch_penalty> <insertion_penalty> <deletion_penalty>"
    echo "Given a sam file, extracts cigar string and calculates mapping score based on score parameters. "
    echo "this script should be in same directory as cigar_to_score.cpp, which will be called by this script automatically"
    echo "match, mismatch, insertion and deletion penalty are hard coded in the cpp file"
    exit
fi

FILENAME=$1
#M_SCORE=$2
#X_SCORE=$3
#I_SCORE=$4
#D_SCORE=$5

NUM_LINES=$(wc -l $FILENAME | awk '{print $1}')
STARTING_LINE=$(awk '{counter +=1 ; if ($1=="@PG") {print counter}}' $FILENAME)
TAIL_LINE_NUM=$(($NUM_LINES - $STARTING_LINE))

tail -n $TAIL_LINE_NUM $FILENAME | awk '{printf "%s %d\n", $6, length($10)}' > cigars.txt &&

g++-8 cigar_to_score.cpp -o cigar_to_score.out && 

./cigar_to_score.out < cigars.txt > scores.txt &&

awk ' {split($0, array, " "); qlen = array[2]; score = array[1]; avgscore=score/qlen; counter += 1; sum += avgscore; sumsqr += avgscore^2} END \
           {printf "Avg score: %f , std: %f out of %d \n", sum/counter, sqrt((sumsqr/counter)-((sum/counter)^2)), counter}' scores.txt
