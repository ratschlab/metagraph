#!/bin/bash

set -e

if [ -z "$6" ]
then
    echo "Usage: $0 <K> <fasta> <uuid> <outbase> <threads> <mem>"
    exit 1
fi

K=$1
fasta=$2
uuid=$3
outbase=$4
threads=$5
mem=$6

## source paths
# metagraph=
# basedir=
. ../paths.sh

if [ ! -f ${outbase}.dbg ]
then 
    ${metagraph} build -v -p $threads -k $K -o ${outbase} --mode canonical --mem-cap-gb $((${mem} / 1250)) --count-kmers ${fasta}
fi

if [ ! -f ${outbase}.clean.fasta.gz ]
then 
    ${metagraph} clean -v --to-fasta -p $threads --prune-tips $(($K * 2)) --prune-unitigs 0 --fallback 3 -o ${outbase}.clean ${outbase}.dbg
fi

rm ${outbase}.dbg* ${outbase}.edgemask

