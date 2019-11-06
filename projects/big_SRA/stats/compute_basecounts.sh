#!/bin/bash

set -e

#module load seqtk
metadata=$(pwd)/../files_WGS_sorted.txt
while IFS=' ' read -r -a line || [[ -n "$line" ]]
do
    fq=${line[0]}
    out=${fq}.fqchk.txt.gz
    if [ ! -f ${out} ]
    then
        echo "module load seqtk; seqtk fqchk -q0 $fq | gzip > $out" | bsub -o /dev/null -M 2000 -We 12:00 -n 1 -J stats
        #echo $fq
        #seqtk fqchk -q0 $fq | gzip > $out
    fi
done < <(cat $metadata)
