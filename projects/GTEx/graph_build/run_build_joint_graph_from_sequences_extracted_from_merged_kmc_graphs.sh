#!/bin/bash

set -e

K=41
threads=24
basedir=/cluster/work/grlab/projects/metagenome/data
kmcdir=${basedir}/gtex/output_k${K}_merged

seqs=""
logfile=${kmcdir}/all_merged_from_seq_k${K}.run.log
rm -f ${logfile}
/usr/bin/time -v ../../metagraph/build/metagengraph build -p $threads --canonical -v -k $K -o ${kmcdir}/all_merged_k${K} ${kmcdir}/*.fasta.gz | tee > ${logfile} 2>&1
