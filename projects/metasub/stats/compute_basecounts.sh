#!/bin/bash

set -e

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv
while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    fqs=${line[40]}
    if [ "$fqs" == "local_raw" ]
    then
        continue
    fi
    for fq in $(echo $fqs | tr ':' ' ')
    do 
        out=${fq}.fqchk.txt.gz
        if [ ! -f ${out} ]
        then
            echo "module load seqtk; seqtk fqchk -q0 $fq | gzip > $out" | bsub -o /dev/null -M 2000 -W 1:00 -n 1 -J stats
        fi
    done
done < <(cat $metadata)
