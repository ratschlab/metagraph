#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 <batch_name>"
    echo "report avg score, std score and runtime from a batch of experiments in separate files easy to copy for plots"
    echo "batch_name should end in either short_reads or long_reads"
    echo "this script should be in same directory as metagraph binary"
    exit
fi

BATCH_NUM=$1
FUNCTION[0]="total_path_score"
FUNCTION[1]="normalized_path_score"
FUNCTION[2]="path_num_matches"

echoerr() { echo "$@" 1>&2; }


rm runtime_tmp.txt
touch runtime_tmp.txt

echo "Extracting runtime"

for FUNC_IND in 0 1 2 
do 
    for USE_CSSW in 1 0
    do
        printf "[" >> runtime_tmp.txt
        for GRAPH_IDX in 0 1 2 3
        do
            printf "[" >> runtime_tmp.txt
            for P in 1 10 100 1000
            do
                EXPERIMENT_DIR="experiments/experiment_batch_${BATCH_NUM}_${FUNCTION[$FUNC_IND]}_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX}"
                if [ ! -d $EXPERIMENT_DIR ]; then
                    echo "Skipping over $EXPERIMENT_DIR. No such experiment."
                    continue
                fi
                if [ ! find $EXPERIMENT_DIR -name "lsf*" > /dev/null 2>&1 ]; then
                    echo "Skipping over $EXPERIMENT_DIR. Experiment in progress."
                    continue
                fi
                LSF_FILE=$(find $EXPERIMENT_DIR -name "lsf*")
                echoerr $EXPERIMENT_DIR
                grep "Exited with exit code" $LSF_FILE && echo "Experiment failed: $EXPERIMENT_DIR"
                grep "CPU time" $LSF_FILE | awk '{line=$0; n=split(line, array, " "); printf "%f, ", array[4];}' >> runtime_tmp.txt
            done
            printf "]\n" >> runtime_tmp.txt
        done
        printf "]\n" >> runtime_tmp.txt
    done
done

cat runtime_tmp.txt | sed ' s/, ]/],/' > runtime_$BATCH_NUM.txt
rm runtime_tmp.txt

echo "Extracting mean scores"

rm scores.txt
touch scores.txt

for FUNC_IND in 0 1 2 
do 
    for USE_CSSW in 1 0 
    do
        printf "[" >> scores.txt
        for GRAPH_IDX in 0 1 2 3
        do
            printf "[" >> scores.txt
            for P in 1 10 100 1000
            do
                EXPERIMENT_DIR="experiments/experiment_batch_${BATCH_NUM}_${FUNCTION[$FUNC_IND]}_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX}"
                if [ ! -d $EXPERIMENT_DIR ]; then
                    echo "Skipping over $EXPERIMENT_DIR. No such experiment."
                    continue
                fi
                echoerr $EXPERIMENT_DIR &&
                ./experiments/alignment_score_tsv.sh $EXPERIMENT_DIR/alignment.sam > score.txt &&
                echoerr $(cat score.txt) &&
                awk '{line=$0; split(line, array, " "); printf "%f, ", array[3]}' score.txt >> scores.txt
                printf "unmapped reads in $EXPERIMENT_DIR: \t"
                grep "$(printf '\t')0" $EXPERIMENT_DIR/alignment.sam | wc -l
                 
            done
            printf "]\n" >> scores.txt
        done
        printf "]\n" >> scores.txt
    done
done

cat scores.txt | sed ' s/, ]/],/' > scores_avg_$BATCH_NUM.txt

echo "Extracting std scores"

rm scores.txt
touch scores.txt

for FUNC_IND in 0 1 2 
do 
    for USE_CSSW in 1 0 
    do
        printf "[" >> scores.txt
        for GRAPH_IDX in 0 1 2 3
        do
            printf "[" >> scores.txt
            for P in 1 10 100 1000
            do
                EXPERIMENT_DIR="experiments/experiment_batch_${BATCH_NUM}_${FUNCTION[$FUNC_IND]}_P_${P}_cssw_${USE_CSSW}_graph_${GRAPH_IDX}"
                if [ ! -d $EXPERIMENT_DIR ]; then
                    echo "Skipping over $EXPERIMENT_DIR. No such experiment."
                    continue
                fi
                echoerr $EXPERIMENT_DIR &&
                ./experiments/alignment_score_tsv.sh $EXPERIMENT_DIR/alignment.sam > score.txt &&
                echoerr $(cat score.txt) &&
                awk '{line=$0; split(line, array, " "); printf "%f, ", array[6]}' score.txt >> scores.txt
            done
            printf "]\n" >> scores.txt
        done
        printf "]\n" >> scores.txt
    done
done

cat scores.txt | sed ' s/, ]/],/; s/]\n]/]]/' > scores_std_$BATCH_NUM.txt
rm score.txt
rm scores.txt
