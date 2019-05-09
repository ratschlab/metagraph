#!/usr/bin/bash

home="/cluster/home/mikhaika"

KMC="$home/projects2014-metagenome/metagraph/build/KMC/kmc"
exe="$home/projects2014-metagenome/metagraph/build/metagengraph_DNA"


if [ $# -ne 1 ]; then
    echo "Usage $0 <fasta.gz>"
    exit 1
fi

FILE="$1"

threshold=0
K=20
num_threads=10

if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

if (( $(stat --printf="%s" $FILE) / (1<<20) < 200 )); then
  cutoff=2
elif (( $(stat --printf="%s" $FILE) / (1<<20) < 500 )); then
  cutoff=3
elif (( $(stat --printf="%s" $FILE) / (1<<20) < 1200 )); then
  cutoff=5
elif (( $(stat --printf="%s" $FILE) / (1<<20) < 2000 )); then
  cutoff=10
else
  cutoff=15
fi

cutoff=$((cutoff-1))
#if [ $cutoff -eq 0 ]; then
#  exit 0
#fi

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

if [ -f $FILE.filter_k$((K-1))_s${cutoff}_${threshold} ] || [ -f $FILE.filter_k$((K-1))_s${cutoff} ]; then
    echo "Already filtered"
else
    mkdir -p "$FILE.cache"
    /usr/bin/time -v $KMC -k$K -m10 -ci1 -fm -t$num_threads $FILE $FILE.kmc $FILE.cache
    rm -r "$FILE.cache"

    /usr/bin/time -v $exe filter -v -p $num_threads -k $((K-1)) --kmc --min-count $((cutoff+1)) --filter-thres $threshold $FILE

    rm $FILE.kmc.*
fi

bsub -J "generate_$(basename $FILE)"_$cutoff -W 48:00 -n 1 -R "rusage[mem=1000]" -o /dev/null \
  "/usr/bin/time -v $exe filter -v --generate-fastq -k $((K-1)) --min-count $((cutoff+1)) --filter-thres $threshold $FILE"

