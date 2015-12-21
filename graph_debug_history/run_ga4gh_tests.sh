#!/bin/bash

set -e

basedir=/cbio/grlab/projects/metagenome/ga4gh
outdir=${basedir}/results

metagraph=/cbio/grlab/home/akahles/git/projects/2014/metagenome/seqan_graph/metagraph
libs=/cbio/grlab/home/akahles/software/libmaus2-lib/lib


for dd in $(find ${basedir}/data/hgvm -maxdepth 1 -mindepth 1 -type d)
do
    if [ "$(basename $dd)" == "CENX" ]
    then
        continue
    fi

    currOutdir=${outdir}/$(basename $dd)/k31
    mkdir -p $currOutdir
    logfile=${currOutdir}/run.log
    echo "cd $currOutdir; export -p LD_LIBRARY_PATH=${libs}:$LD_LIBRARY_PATH; $metagraph -k 31 -v ${dd}/ref.fa ${dd}/[OG]*.fa" | qsub -l nodes=1:ppn=1,walltime=12:00:00,vmem=4G,pmem=4G,mem=4G -N metagraph -j oe -o $logfile

    currOutdir=${outdir}/$(basename $dd)/k63
    mkdir -p $currOutdir
    logfile=${currOutdir}/run.log
    echo "cd $currOutdir; export -p LD_LIBRARY_PATH=${libs}:$LD_LIBRARY_PATH; $metagraph -k 63 -v ${dd}/ref.fa ${dd}/[OG]*.fa" | qsub -l nodes=1:ppn=1,walltime=12:00:00,vmem=4G,pmem=4G,mem=4G -N metagraph -j oe -o $logfile
done
