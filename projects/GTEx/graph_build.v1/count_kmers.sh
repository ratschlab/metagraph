#!/usr/bin/bash


KMC="$(dirname ${BASH_SOURCE[0]})/../../../metagraph/build/KMC/kmc"


if [ $# -ne 4 ]; then
    echo -e "Usage:\n$0 <k> <file_with_files> <outdir> <sampleID>" >&2
    exit 1
fi

K="$1"
FILES="$2"
OUT="$3"
SAMPLE="$4"
cutoff=2
num_threads=4

if [ ! -f $FILES ]; then
    echo "File" $FILES "does not exist" >&2
    exit 1
fi

if [ $cutoff -eq 0 ]; then
  echo "Error: cutoff is too small" >&2
  exit 0
fi

#filename="$(basename $FILE)"
mkdir -p "/scratch/$SAMPLE.cache"
/usr/bin/time -v $KMC -k$K -m10 -ci$cutoff -fq -t$num_threads "@"$FILES $OUT/${SAMPLE}.k${K} /scratch/$SAMPLE.cache
rm -r "/scratch/$SAMPLE.cache"
