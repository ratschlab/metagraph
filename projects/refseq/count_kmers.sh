#!/usr/bin/bash


KMC="$(dirname ${BASH_SOURCE[0]})/../../metagraph/build/KMC/kmc"


if [ $# -ne 3 ]; then
    echo -e "Usage:\n$0 <k> <file.fasta.gz> <outdir>" >&2
    exit 1
fi

K="$1"
FILE="$2"
OUT="$3"
cutoff=1
num_threads=10

if [ ! -f $FILE ]; then
    echo "File" $FILE "does not exist" >&2
    exit 1
fi

if [ $cutoff -eq 0 ]; then
  echo "Error: cutoff is too small" >&2
  exit 0
fi

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected" >&2
    exit 1
fi

filename="$(basename $FILE)"
mkdir -p "/scratch/$filename.cache"
/usr/bin/time -v $KMC -k$K -m10 -ci$cutoff -b -fm -t$num_threads $FILE $OUT/$filename /scratch/$filename.cache
rm -r "/scratch/$filename.cache"
