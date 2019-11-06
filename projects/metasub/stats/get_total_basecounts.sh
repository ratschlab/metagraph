#!/bin/bash

set -e

fastqdir=/cluster/work/grlab/projects/metagenome/raw_data/metasub/wasabi_raw_only
outfile=${fastqdir}.cnt.stats

if [ ! -f ${outfile} ]
then
    find ${fastqdir} -name \*.fqchk.txt.gz -exec zcat {} \; | grep ALL | cut -f 2 >  $outfile
fi

### display total count
awk '{sum += $1} END {print sum}' ${outfile}
