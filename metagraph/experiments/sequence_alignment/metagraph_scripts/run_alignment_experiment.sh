#!/bin/bash

if [ $# -lt 9 ];
then
    echo "Usage: run_alignment_experiment.sh <experiment_name> <data_file> <graph_file> <graph_type> <queue_size> <path_comparison_function> <use cssw lib(set 0 or 1)> <num_cores> <memory per core> <extra_flags>"
    echo "this script should be in same directory as metagraph binary"
    exit
fi
EXPERIMENT_NAME=$1
DIR_NAME=experiments/$1
DATA_FILE=$2
GRAPH=$3
GRAPH_TYPE=$4
NUM_PATHS=$5
PATH_COMPARE_FUNC=$6
USE_CSSW=$7
NUM_CORES=$8
MEMORY=$9

SCORE_SCRIPT="./experiments/alignment_score_tsv.sh"

EXTRA_FLAGS="${10} --discard-similar-paths"
if [ $USE_CSSW -ne 0 ]; then
    EXTRA_FLAGS="$EXTRA_FLAGS --align-using-cssw-library"
fi

if [ -d ${DIR_NAME} ];
then
    echo "Experiment already exists"
    exit
fi
mkdir $DIR_NAME
pushd $DIR_NAME

git log -n 5 > git_log
git diff > git_diff

bsub -J "align" -R"rusage[mem="$MEMORY"]" -R"select[hname='lm-fat-001']" -q normal.24h -n $NUM_CORES /cluster/home/jsara/projects2014-metagenome/metagraph/build_profile/metagengraph align -i $GRAPH $DATA_FILE --graph $GRAPH_TYPE --align-num-paths $NUM_PATHS -o alignment.sam --align-path-comparison-function $PATH_COMPARE_FUNC $EXTRA_FLAGS
#bsub -J "align" -R"rusage[mem="$MEMORY"]" -R"select[hname='lm-s4-051']" -q normal.24h -n $NUM_CORES /cluster/home/jsara/projects2014-metagenome/metagraph/build_profile/metagengraph align -i $GRAPH $DATA_FILE --graph $GRAPH_TYPE --align-num-paths $NUM_PATHS -o alignment.sam --align-path-comparison-function $PATH_COMPARE_FUNC $EXTRA_FLAGS
echo "Results should be stored in ${DIR_NAME}."
popd
