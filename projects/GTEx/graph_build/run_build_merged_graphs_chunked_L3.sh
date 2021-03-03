#!/bin/bash

set -e

## source paths
# set metagraph
# set basedir
. ../paths.sh

K=41

mem=250000
threads=6
threads=18
pmem=$(($mem / $threads))
seqdir=${basedir}/graphs/gtex/output_k${K}_trimmed_clean
outdir=${basedir}/graphs/gtex/output_k${K}_trimmed_clean_graph_chunked
mkdir -p $outdir

if [ ! -f ${seqdir}_all_files.txt ]
then
    find ${seqdir} -name \*.fasta.gz > ${seqdir}_all_files.txt
fi

#for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}
for F in \\\$\\\$\\\$
do
    echo "cat ${seqdir}_all_files.txt | /usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_merged_k${K} --suffix $F 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J gt_g_${F} -oo ${outdir}/build_k${K}_$F.lsf.log -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done

