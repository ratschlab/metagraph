#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/results/gtex
all_files=${basedir}/metadata/SraRunTable_20180218.txt

for tissue in $(cut -f 24 $all_files | sort -u | tr ' ' '_')
do
    if [ "$tissue" == "histological_type" ]
    then
        continue
    fi
    grep -e "$(echo $tissue | tr '_' ' ')" $all_files > ${basedir}/metadata/by_tissue/SraRunTable_20180218.${tissue}.txt
done
