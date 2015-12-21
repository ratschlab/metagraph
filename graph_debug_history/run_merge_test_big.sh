#!/bin/bash

set -e

basedir=/cbio/grlab/projects/metagenome/bacteria/results/27
basedir2=/cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/

./metagraph -v -k 27 -O test_merge2_big_out ${basedir2}/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377/*.fna ${basedir2}/Acetobacter_pasteurianus_386B_uid214433/*.fna 
md5sum test_merge2_big_out.W.dbg
md5sum test_merge2_big_out.l.dbg

#./metagraph -v -m ${basedir}/Acetobacter_pasteurianus_IFO_3283_01_42C_uid158377,${basedir}/Acetobacter_pasteurianus_386B_uid214433 -O test_merge_big_out BBB
md5sum test_merge_big_out.W.dbg
md5sum test_merge_big_out.l.dbg
