#!/bin/bash

set -e

k=$1

if [ -z "$k" ]
then
    echo "Usage: $0 k" 
    exit 1
fi

total=1000

for i in $(seq $(($total + 1)) $(($total + 1)))
do

    head -n 1 /cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377/NC_017104.fna > medium1_test.fa
    head -n $i /cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377/NC_017104.fna | tail -n $total >> medium1_test.fa
    head -n 1 /cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377/NC_017105.fna > medium2_test.fa
    head -n $i /cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377/NC_017105.fna | tail -n $total >> medium2_test.fa

    #./metagraph -v -k $k -O medium1_test medium1_test.fa
    #./metagraph -v -k $k -O medium2_test medium2_test.fa

    time ./metagraph -v -m medium1_test,medium2_test -O test_merge_medium_out BBB
    #md5sum test_merge_medium_out.W.dbg | cut -f 1 -d ' ' > tut1 
    #md5sum test_merge_medium_out.l.dbg

    ./metagraph -v -k $k -O test_merge2_medium_out medium1_test.fa medium2_test.fa 
    #md5sum test_merge2_medium_out.W.dbg | cut -f 1 -d ' ' > tut2
    #md5sum test_merge2_medium_out.l.dbg

    ./metagraph -v -c test_merge_medium_out,test_merge2_medium_out BBB
done

