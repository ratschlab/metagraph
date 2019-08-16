#!/usr/bin/bash

KMC="$(dirname ${BASH_SOURCE[0]})/../build/KMC/kmc"
exe="$(dirname ${BASH_SOURCE[0]})/../build/metagraph"


if [ $# -ne 3 ]; then
    echo -e "Usage:\n$0 <fastq.gz> <cutoff> <threshold>" >&2
    exit 1
fi

FILE="$1"
cutoff=$2
threshold=$3
K=21
num_threads=10

if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

if [ -f $FILE.filter_k$((K-1))_s${cutoff}_${threshold} ]; then
    echo "Already filtered"
    exit 0
fi

mkdir -p "$FILE.cache"
/usr/bin/time -v $KMC -k$K -m10 -ci1 -fq -t$num_threads $FILE $FILE.kmc $FILE.cache
rm -r "$FILE.cache"

/usr/bin/time -v $exe filter -v -p $num_threads -k $((K-1)) --kmc --min-count $((cutoff+1)) --filter-thres $threshold $FILE

rm $FILE.kmc.*
