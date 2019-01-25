#!/usr/bin/bash

exe="$(dirname ${BASH_SOURCE[0]})/../build/metagengraph"


if [ $# -ne 2 ] || [[ ${1: -8} != '.kmc_suf' ]]; then
    echo -e "Usage:\n$0 <genomic.kmc_suf> <k>" >&2
    exit 1
fi

if [ ! -f $1 ]; then
    echo "File does not exist"
    exit 1
fi

FILE="${1%.kmc_suf}"
K=$2
num_threads=10

/usr/bin/time -v $exe build -v -p $num_threads -k $K -o $FILE --kmc $FILE
/usr/bin/time -v $exe transform -v --to-fasta --contigs -p $num_threads -o ${FILE}.contigs $FILE
