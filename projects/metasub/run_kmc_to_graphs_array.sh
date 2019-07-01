#!/bin/bash

set -e

K=19

mem=50000
threads=4
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metasub/graphs/output_k${K}
mkdir -p $outdir

filelist=${outdir}/metasub_kmc_k${K}.txt
if [ ! -f ${filelist} ]
then
    ls -1 ${basedir}/kmc_counts/output_k${K}/*.kmc_suf > ${filelist}
fi

echo "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
      if [ -f ${outdir}/\${x%.kmc_suf}.dbg ]; then exit 0; fi; \
      ../../metagraph/build/metagraph build -p $threads -k $K -o ${outdir}/\${x%.kmc_suf} --kmc \${kmc_file}" | bsub -J contigs_k${K}[1-4173]%400 -o ${basedir}/graph_build_${K}.lsf -We 20:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]"
