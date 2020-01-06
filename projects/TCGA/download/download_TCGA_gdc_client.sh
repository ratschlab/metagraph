#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

basedir=/cluster/work/grlab/projects/metagenome/raw_data/tcga
datadir=${basedir}/data/${tissue}
gdc=/cluster/home/akahles/software/gdc-client/gdc-client
keys=/cluster/home/akahles/.keys
mkdir -p ${datadir}

cnt=0
cnt_done=0
while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    aid=${line[0]}
    if [ "$aid" == "id" ]
    then 
        continue
    fi

    if [ ! -d ${datadir}/${aid} ]
    then
        echo "Downloading to $datadir $aid"
        $gdc download -t ${keys}/gdc-user-token.2020-03-20T16_36_32.869Z.txt -n 2 -d $datadir $aid & #--http-chunk-size 65536 $uuid
        cnt=$(($cnt + 1))
    else
        echo " ${datadir}/${aid} already exists"
        cnt_done=$(($cnt_done + 1))
    fi
    if [ "$cnt" == "8" ]
    then
        wait
        cnt=0
    fi
done < <(cat ${basedir}/metadata/gdc_manifest.2020-02-22.${tissue}.txt | tr $'\t' ',')
echo $cnt_done samples were already finished
