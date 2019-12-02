#!/bin/bash

set -e

fastqdir=/cluster/work/grlab/projects/metagenome/data/MetaGut/nobackup/human_gut_sra
outfile=${fastqdir}.WGS.cnt.stats

if [ ! -f ${outfile} ]
then
    find ${fastqdir} -name \*.fqchk.txt.gz -exec zcat {} \; | grep ALL | cut -f 2 >  $outfile
fi

### display total count
awk '{sum += $1} END {print sum}' ${outfile}
