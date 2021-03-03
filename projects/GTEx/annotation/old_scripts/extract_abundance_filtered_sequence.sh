#!/bin/bash

set -e

### source paths
. ../paths.sh

if [ -z "$3" ]
then
    echo "Usage: $0 <kmc_file> <out_dir> <K>"
    exit 1
fi
kmc_file=$1
out_dir=$2
K=$3

outbase=$(basename $kmc_file)
outbase=${outbase%.kmc_suf}

tmp=${out_dir}/tmp.${outbase}.${RANDOM}
mkdir -p $tmp

quantiles=( "0.0" "0.2" "0.4" "0.6" "0.8" "1.0" )

for i in $(seq 0 $((${#quantiles[@]} - 2)))
do
    # already complete
    if [ -f "${out_dir}/${outbase}.exp_level${i}.sequences.fasta.gz" ]
    then
        continue
    fi

    echo processing level $i: ${quantiles[$i]} to ${quantiles[$(($i + 1))]}
    out_filt_base=${tmp}/filtered.level${i}

    # generate dbg
    echo building graph
    $metagraph build -p 2 -k $K -o ${out_filt_base} --mode canonical --min-count-q ${quantiles[$i]} --max-count-q ${quantiles[$(($i + 1))]} ${kmc_file}

    # extract sequence from graph
    echo extracting sequencinges:
    $metagraph transform --to-fasta -o ${out_dir}/${outbase}.exp_level${i}.sequences ${out_filt_base}.dbg

    # clean up
    rm ${tmp}/*
done
rm -rf $tmp
