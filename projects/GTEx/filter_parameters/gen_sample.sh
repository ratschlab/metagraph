#!/bin/bash

set -e

out=samples_parameterstudy.tsv
rm -rf $out
lines="192 468 1021 1408 2843 4999 5188 6112 8038 9300"
for i in ${lines} #$(seq 1 10)
do
    grep -w RNA-Seq /cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt | sed -n ${i}p >> $out
done
