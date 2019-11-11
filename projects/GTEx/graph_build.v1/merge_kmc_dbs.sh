#!/bin/bash

set -e 

if [ -z "$2" ]
then
    echo "Usage: $0 <sample_count> <threshold> [<FROM>]"
    exit 1
fi
TOP="$1"
THRESH="$2"
FROM=1
if [ ! -z "$3" ]
then
    FROM="$3"
fi
TO=$((${FROM} + ${TOP}))

threads=8
K=41

module load kmc

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/
datadir=${basedir}/output_k${K}
outdir=${basedir}/output_k${K}_merged
mkdir -p $outdir

cd $outdir

### generate input file
out=kmc_merge_file.${FROM}_${TO}.t${THRESH}.txt
echo "INPUT:" > $out
ops=""
for fname in $(ls -1 ${datadir}/*.kmc_suf | tail -n+${FROM} | head -n $TOP)
do
    fbase=$(basename $fname | cut -f 1 -d '.')
    echo "${fbase} = ${fname%.kmc_suf}" >> $out
    ops="${ops} + ${fbase}"
done
echo "OUTPUT:" >> $out
echo "merged.${FROM}_${TO}.t${THRESH} = ${ops# +}" >> $out
echo "OUTPUT_PARAMS:" >> $out
echo "-ci${THRESH}" >> $out
echo "-cs10" >> $out

### run kmc_tools
kmc_tools -t${threads} complex $out
