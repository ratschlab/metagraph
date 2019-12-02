#!/bin/bash

set -e

## source paths
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

#mem=60000
#mem=100000
#mem=240000
mem=360000
mem=420000
mem=600000
threads=6
threads=12
pmem=$((${mem} / ${threads}))

outdir=${basedir}/metagut/graphs/output_k${K}
logdir=${basedir}/metagut/graphs/output_k${K}_logs
mkdir -p $outdir
mkdir -p $logdir

#metadata=$(pwd)/../files_WGS_sorted.2000.txt
#metadata=$(pwd)/../files_WGS_sorted.4000.txt
#metadata=$(pwd)/../files_WGS_sorted.8000.txt
#metadata=$(pwd)/../files_WGS_sorted.12000.txt
#metadata=$(pwd)/../files_WGS_sorted.20000.txt
metadata=$(pwd)/../files_WGS_sorted.txt


while IFS=' ' read -r -a line || [[ -n "$line" ]]
do
    fasta=${line[0]}
    uuid=$(basename $fasta | cut -f 1 -d '.')

    outbase=${outdir}/${uuid}
    logfile1=${logdir}/graph_build.${uuid}.lsf.log
    logfile2=${logdir}/graph_build.${uuid}.log

    if [ -f ${outbase}.clean.fasta.gz ]
    then
        echo ${outbase} complete
        continue
    fi
    
    echo "cd $(pwd); /usr/bin/time -v ./build_clean_graphs.sh $K $fasta $uuid $outbase $threads $mem | tee $logfile2" | bsub -g /akahles/metagut -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J metagut_k${K} -oo $logfile1
done < <(cat $metadata)
