#!/bin/bash

set -e

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

module load kmc
metagraph=../../metagraph/build/metagengraph

lower_bounds=(2 3 5 17 255)
upper_bounds=(2 4 16 254 10000)

for i in $(seq 0 $((${#lower_bounds[@]} - 1)))
do
    # already complete
    if [ -f "${out_dir}/${outbase}.exp_level${i}.sequences.fasta.gz" ]
    then
        continue
    fi

    # generate filtered KMC dump
    echo processing level $i: ${lower_bounds[$i]} to ${upper_bounds[$i]}
    out_filt_base=${tmp}/kmc_filt.level${i}
    kmc_tools transform ${kmc_file%.kmc_suf} reduce $out_filt_base -ci${lower_bounds[$i]} -cx${upper_bounds[$i]}

    # generate dbg from filtered KMC
    echo building graph
    $metagraph build -p 2 -k $K -o ${out_filt_base} --kmc ${out_filt_base}.kmc_suf

    # extract sequence from graph
    echo extracting sequencinges:
    $metagraph transform --to-fasta -o ${out_dir}/${outbase}.exp_level${i}.sequences ${out_filt_base}.dbg

    # clean up
    rm ${tmp}/*
done
rm -rf $tmp
