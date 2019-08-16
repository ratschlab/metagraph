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
mem=32000

if [ ! -f ${FILE}.k${K}.dbg ]
then
    /usr/bin/time -v $exe build -v -p $num_threads -k $K -o ${FILE}.k$K --kmc $FILE || exit 1
fi
if [ ! -f ${FILE}.k$K.contigs.fa.gz ]
then
    echo "/usr/bin/time -v $exe assemble -v --unitigs -p 1 -o ${FILE}.k$K.contigs ${FILE}.k$K" | bsub -M $mem -n 1 -We 12:00 -R "rusage[mem=${mem}]" -J kmc_msub2 -o ${FILE}.kmc2contS2.k${K}.cluster.log
fi
