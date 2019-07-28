#!/bin/bash

set -e

K=19

outdir=/cluster/work/grlab/projects/metagenome/data/metasub/column_annotation_k${K}
graph=${outdir}/metasub_wasabi_complete.k${K}.bitmapdbg
mem=50000
mkdir -p ${outdir}

if [ ! -f ${graph} ]
then
    echo building graph for k $K
    $(pwd)/../../metagraph/build/metagraph build --complete --graph bitmap --canonical -k $K -o ${graph%.bitmapdbg}
fi

for fname in /cluster/work/grlab/projects/metagenome/data/metasub/contigs/metasub_k${K}/*.fasta.gz
do
    sample=$(basename $fname | cut -f 1 -d '.')
    outbase=${outdir}/metasub_wasabi_complete.k${K}.f2.${sample}
    logfile=${outbase}.cluster.log
    if [ ! -f ${outbase}.column.annodbg ]
    then
        #echo "/usr/bin/time -v $(pwd)/../../metagraph/build/metagraph annotate -i $graph -o $outbase --anno-filename $fname > $logfile 2>&1" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 4:00 -J anno -o /dev/null
        echo "/usr/bin/time -v $(pwd)/../../metagraph/build/metagraph annotate -i $graph -o $outbase --anno-filename $fname" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 4:00 -J anno -o $logfile
    else
        echo "$sample already processed"
    fi
done
