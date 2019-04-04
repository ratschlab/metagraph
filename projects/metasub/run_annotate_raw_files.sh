#!/bin/bash

set -e

K=15

outdir=/cluster/work/grlab/projects/metagenome/data/metasub/column_annotation_raw_k${K}
graph=${outdir}/metasub_wasabi_complete.k${K}.bitmapdbg
mem=50000
mkdir -p $outdir

for fname in $(find /cluster/work/grlab/projects/metagenome/raw_data/metasub/wasabi_raw_only -name \*_1.fastq.gz)
do
    sample=$(basename $fname | cut -f 1 -d '.' | sed -e "s/_1$//g")
    outbase=${outdir}/metasub_wasabi_complete.k${K}.f1.${sample}
    logfile=${outbase}.cluster.log
    fname2=${fname%_1.fastq.gz}_2.fastq.gz
    if [ ! -f ${outbase}.column.annodbg ]
    then
        echo "/usr/bin/time -v $(pwd)/../../metagraph/build/metagengraph annotate -i $graph -o $outbase --anno-filename $fname $name2 > $logfile 2>&1" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 4:00 -J anno -o /dev/null
    else
        echo "$sample already processed"
    fi
done
