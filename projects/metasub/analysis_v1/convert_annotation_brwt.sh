#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data
annograph=/cluster/home/akahles/git/software/genome_graph_annotation/build/annograph
K=17
F=2
mem=50000
threads=30

echo "/usr/bin/time -v $annograph transform_anno -v -o ${basedir}/metasub_wasabi_graph_anno.k${K}.f${F} --anno-type brwt --greedy ${basedir}/metasub_wasabi_graph_anno.k${K}.f${F}.column.annodbg -p 60 2>&1" | bsub -J convert_metasub_to_brwt_pm -oo conversion_to_brwt.f${F}.k${K}.lsf -W 150:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]" 
