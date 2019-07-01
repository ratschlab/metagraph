#!/bin/bash

set -e

K=41
threads=16
basedir=/cluster/work/grlab/projects/metagenome/data
kmcdir=${basedir}/gtex/output_k${K}_merged

rm -f ${kmcdir}/all_merged_k${K}.run.log
for fname in $(ls -1 ${kmcdir}/*.kmc_suf)
do
    fbase=$(basename $fname)
    if [ ! -f ${kmcdir}/${fbase%.kmc_suf}.dbg ]
    then
        ../../../metagraph/build/metagengraph build -p $threads -v -k $K -o ${kmcdir}/${fbase%.kmc_suf} --kmc ${fname} >> ${kmcdir}/all_merged_k${K}.run.log
    fi
done
