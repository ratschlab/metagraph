#!/bin/bash

set -e

K=41
threads=16
basedir=/cluster/work/grlab/projects/metagenome/data
kmcdir=${basedir}/gtex/output_k${K}_merged
mem=50000
pmem=50000
threads=1

logfile=${kmcdir}/all_merged_extract_seqs_k${K}.run.log
rm -f ${logfile}
for fname in $(ls -1 ${kmcdir}/merged*.dbg)
do
    fbase=$(basename $fname)
    if [ ! -f ${kmcdir}/${fbase%.dbg}.fasta.gz ]
    then
        echo "$(pwd)/../../../metagraph/build/metagengraph transform --to-fasta -v -o ${kmcdir}/${fbase%.dbg}.sequences ${fname}" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -n 1 -J gtex_seq -o /dev/null
    fi
done
