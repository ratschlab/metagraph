#!/bin/bash

set -e

## source paths
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$2" ]
then
    echo "Usage: $0 <K <tissue>>"
    exit 1
fi
K=$1
tissue=$2

#mem=125000
mem=250000
#mem=350000
#mem=500000
threads=4
pmem=$((${mem} / ${threads}))

countdir=${basedir}/tcga/output_k${K}_trimmed/${tissue}
outdir=${basedir}/tcga/output_k${K}_trimmed_clean/${tissue}
mkdir -p $outdir

for kmc_file in ${countdir}/*.kmc_suf
do
    x=$(basename $kmc_file)
    outbase=${outdir}/${x%.kmc_suf}
    if [ -f ${outbase}.clean.0.000000.0.200000.fasta.gz ]
   # if [ -f ${outbase}.clean.fasta.gz ]
    then
        echo $x complete
        continue
    fi
    #echo "${metagraph} build -p $threads -v -k $K -o ${outbase} --mode canonical --count-kmers ${kmc_file}; ${metagraph} clean --to-fasta -p ${threads} --prune-tips $(($K * 2)) --prune-unitigs 0 --fallback 3 --unitigs --count-bins-q '0 0.2 0.4 0.6 0.8 1.0' -o ${outbase}.clean ${outbase}.dbg; rm ${outbase}.dbg* ${outbase}.edgemask" | bsub -G ms_raets -J tcga_clean_k${K} -oo ${outbase}.lsf.log -We 20:00 -n $threads -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -M ${mem}
    echo "${metagraph} build -p $threads -v -k $K -o ${outbase} --mode canonical --count-kmers ${kmc_file}; ${metagraph} clean --to-fasta -p ${threads} --prune-unitigs 0 --fallback 3 --unitigs --count-bins-q '0 0.2 0.4 0.6 0.8 1.0' -o ${outbase}.clean ${outbase}.dbg; rm ${outbase}.dbg* ${outbase}.edgemask" | bsub -G ms_raets -J tcga_clean_k${K} -oo ${outbase}.lsf.log -We 20:00 -n $threads -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -M ${mem}
done
