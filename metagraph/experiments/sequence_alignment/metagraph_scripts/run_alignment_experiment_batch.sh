#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: ./run_alignment_experiment_batch.sh <batch_number> <input_data_file_for_short_reads> <input_data_file_for_long_reads>"
    echo "this script should be in same directory as metagraph binary"
    exit
fi

BATCH_NUM=$1
INPUT_DATA_SHORT=$2
INPUT_DATA_LONG=$3

GRAPHS[0]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/w_variations/hg19_hs37d5_succinct_k_12_w_vcf_all.dbg"
GRAPHS[1]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/w_variations/hg19_hs37d5_succinct_k_15_w_vcf_all.dbg"
GRAPHS[2]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/w_variations/hg19_hs37d5_succinct_k_21_w_vcf_all.dbg"
GRAPHS[3]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/w_variations/hg19_hs37d5_succinct_k_27_w_vcf_all.dbg"

#GRAPHS[0]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/hg19_hs37d5_succinct_k_12.dbg"
#GRAPHS[1]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/hg19_hs37d5_succinct_k_15.dbg"
#GRAPHS[2]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/hg19_hs37d5_succinct_k_21.dbg"
#GRAPHS[3]="/cluster/work/grlab/projects/metagenome/data/alignment/graphs/hg19_hs37d5_succinct_k_27.dbg"


for GRAPH_IDX in 0 1 2 3
do
    GRAPH=${GRAPHS[$GRAPH_IDX]}
    for P in 1 10 100 1000
    do 
        for USE_CSSW in 0 1
        do
            echo "Running experiment short_reads_total_path_score_cssw and num_paths = $P and flag for using cssw: $USE_CSSW" &&
            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_short_reads_total_path_score_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_SHORT $GRAPH succinct $P 0 $USE_CSSW 1 10000
            echo "Running experiment short_reads_normalized_path_score_cssw and num_paths = $P and flag for using cssw: $USE_CSSW" &&
            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_short_reads_normalized_path_score_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_SHORT $GRAPH succinct $P 1 $USE_CSSW 1 10000
            echo "Running experiment short_reads_path_num_matches_cssw and num_paths = $P and flag for using cssw: $USE_CSSW" &&
            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_short_reads_path_num_matches_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_SHORT $GRAPH succinct $P 2 $USE_CSSW 1 10000
            echo "Submitted short reads experiments"
        done
            echo "Running experiment long_reads_total_path_score_cssw and num_paths = $P and flag for using cssw: 1" &&
            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_long_reads_total_path_score_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_LONG $GRAPH succinct $P 0 1 1 10000 --disable-cssw-speedup 
#            echo "Running experiment long_reads_normalized_path_score_cssw and num_paths = $P and flag for using cssw: $USE_CSSW" &&
#            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_long_reads_normalized_path_score_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_LONG $GRAPH succinct $P 1 $USE_CSSW 1 10000 &&
#            echo "Running experiment long_reads_path_num_matches_cssw  and num_paths = $P and flag for using cssw: $USE_CSSW" &&
#            ./run_alignment_experiment.sh experiment_batch_"$BATCH_NUM"_long_reads_path_num_matches_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX} $INPUT_DATA_LONG $GRAPH succinct $P 2 $USE_CSSW 1 10000
    done
done 
echo "All experiments are running :-)"
