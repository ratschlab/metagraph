#!/usr/bin/bash

home="/cluster/home/mikhaika"

KMC="$home/projects2014-metagenome/metagraph/build/KMC/kmc"
exe="$home/projects2014-metagenome/metagraph/build/metagraph_DNA"


if [ $# -ne 3 ]; then
    echo "Usage $0 <fasta.gz> <cutoff> <threshold>"
    exit 1
fi

FILE="$1"
cutoff=$2
threshold=$3
K=25
num_threads=10

if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

cutoff=15

cutoff=$((cutoff-1))

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

if [ -f $FILE.filter_k$((K-1))_s${cutoff}_${threshold} ]; then
    echo "Already filtered"
    exit 0
fi

if [ -f $FILE.filter_k$((K-1))_s${cutoff} ]; then
    echo "Already filtered"
    exit 0
fi

mkdir -p "$FILE.cache"
/usr/bin/time -v $KMC -k$K -m21 -ci1 -fm -t$num_threads $FILE $FILE.kmc $FILE.cache
rm -r "$FILE.cache"

/usr/bin/time -v $exe filter -v -p $num_threads -k $((K-1)) --kmc --min-count $((cutoff+1)) --filter-thres $threshold $FILE

rm $FILE.kmc.*

#bsub -J "generate_$(basename $FILE)"_$cutoff -W 48:00 -n 1 -R "rusage[mem=1000]" -o /dev/null \
#  "/usr/bin/time -v $exe filter -v --generate-fastq -k $((K-1)) --min-count $((cutoff+1)) --filter-thres $threshold $FILE"

