#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 corr file"
    exit
fi

INPUT_FILE=$1

awk '{line = $0; split(line, words, " "); split(words[3], is_mapped_tail, "("); split(is_mapped_tail[2], is_mapped_head, ")"); is_mapped = is_mapped_head[1]; \
    if (is_mapped == 1 ) {split(words[4], score_tail, "("); split(score_tail[2], score_head, ")"); n=split(score_head[1], score, "-"); if(n == 1) {print(score[1])} else {print(score[2])};}}' $INPUT_FILE \
    | awk ' {line = $0; line = line * 2; qlen = 100.0; num=line/qlen; counter += 1; sum += num; sumsqr += num^2} END \
           {printf "Avg score: %f , std: %f out of %d \n", sum/counter, sqrt((sumsqr/counter)-((sum/counter)^2)), counter}'
