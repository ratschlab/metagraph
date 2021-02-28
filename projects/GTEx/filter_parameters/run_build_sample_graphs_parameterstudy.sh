#!/bin/bash

set -e

K=41

mem=250000
threads=4
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data
#outdir=${basedir}/gtex/parameterstudy/output_k${K}
outdir=${basedir}/gtex/parameterstudy_largest/output_k${K}
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current2/metagraph/build/metagraph
#metadata="samples_parameterstudy.tsv"
metadata="samples_parameterstudy.largest.tsv"
kmcdir=${basedir}/gtex/output_k${K}/

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[16]}
    kmc=${kmcdir}/${uuid}.k${K}.kmc_suf

    outfname=${outdir}/${uuid}.dbg
    logfname=${outdir}/${uuid}.lsf.log
    if [ -f ${outfname} ]
    then
        echo "$uuid complete"
        continue
    fi

    echo "/usr/bin/time -v $metagraph build -v -p $threads -k $K -o ${outfname%.dbg} --mode canonical --count-kmers $kmc" | bsub -J mg_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done < <(cat ${metadata} | tr $'\t' ',')
