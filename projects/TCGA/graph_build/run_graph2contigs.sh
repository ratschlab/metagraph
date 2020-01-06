#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

K=31

mem=50000
threads=1
pmem=$((${mem} / ${threads}))

basedir=/cluster/work/grlab/projects/metagenome/data
graphbase=${basedir}/tcga/output_k${K}_trimmed_clean_graph_chunked/${tissue}/graph_merged_k${K}
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

if [ -f ${graphbase}.fasta.gz ]
then
    echo $x complete
    continue
fi
echo "${metagraph} transform --to-fasta --unitigs -o ${graphbase} ${graphbase}.dbg" | bsub -J tc_gr2fa_k${K} -oo ${graphbase}.lsf.log -We 20:00 -n $threads -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -M ${mem}
