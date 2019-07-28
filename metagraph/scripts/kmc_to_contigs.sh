#!/usr/bin/bash

exe="$(dirname ${BASH_SOURCE[0]})/../build/metagraph"


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
num_threads=4

/usr/bin/time -v $exe build -v -p $num_threads -k $K -o $FILE --kmc $FILE
/usr/bin/time -v $exe assemble -v --unitigs -p $num_threads -o ${FILE}.k$K.contigs $FILE
