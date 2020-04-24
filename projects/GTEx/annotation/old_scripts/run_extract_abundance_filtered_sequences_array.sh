#!/bin/bash

### source paths
# basedir=
. ../paths.sh

set -e

mem=30000
mem=80000
mem=150000
threads=2
pmem=$(($mem / $threads))
K=41
outdir=${basedir}/gtex/graphs/output_k${K}_sequences_per_abundance_quantiles
mkdir -p $outdir
logdir=${basedir}/gtex/graphs/output_k${K}_sequences_per_abundance_quantiles_logs
mkdir -p $logdir
filelist=${basedir}/gtex/graphs/gtex_kmc_k${K}.txt
N=9759

echo  "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
$(pwd)/extract_abundance_filtered_sequence.sh \
\${kmc_file} ${outdir} $K" | bsub -J contigs_k${K}[1-${N}]%400 -M $mem -o $logdir/seq_ext_abundances_k${K}_quantiles.%I.lsf -We 20:00 -n $threads -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
