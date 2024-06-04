```bash
DIR=~/metagenome/data/gtex_counts/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
KMC=~/projects/projects2014-metagenome/metagraph/build_release/KMC/kmc;

mkdir $DIR;
mkdir $DIR/logs;
mkdir $DIR/graphs;
mkdir $DIR/unitigs;

find /cluster/work/grlab_shared/GTEx/extract/dbGaP-9608/fastq -name "*_1.fastq.gz" > $DIR/input_fastq.txt

sbatch -J "build_graphs" \
     --array=1-9858 \
     -o $DIR/logs/build_graphs.slog \
     -t 00-7 \
     --cpus-per-task 4 \
     --mem-per-cpu=15G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="file=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${DIR}/input_fastq.txt); \
        file=\${file%_1.fastq.gz}; \
        name=\$(basename \${file}); \
        [ -f $DIR/unitigs/\$name.fasta.gz ] && exit 0; \
        echo \"\${file}_1.fastq.gz\n\${file}_2.fastq.gz\" > $DIR/graphs/\${name}.input; \
        mkdir ~/metagenome/scratch/nobackup/stripe_1/\${name}.kmc_cache; \
        /usr/bin/time -v $KMC -k31 -m40 -sm -ci1 -cs65535 -fq -t4 \
            @$DIR/graphs/\${name}.input \
            $DIR/graphs/\$name \
            ~/metagenome/scratch/nobackup/stripe_1/\${name}.kmc_cache \
        > $DIR/graphs/\$name.kmc.log 2>&1; \
        rm -r ~/metagenome/scratch/nobackup/stripe_1/\${name}.kmc_cache; \
        /usr/bin/time -v $METAGRAPH build -v -k 31 \
            --mode canonical \
            --count-kmers --count-width 16 \
            --mem-cap-gb 30 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 4 \
            -o $DIR/graphs/\$name \
            $DIR/graphs/\$name.kmc_suf > $DIR/graphs/\$name.build.log 2>&1 \
        && rm $DIR/graphs/\$name.kmc_* \
        && /usr/bin/time -v $METAGRAPH clean \
            --to-fasta --primary-kmers \
            --prune-unitigs 0 \
            --fallback 2 \
            -p 4 \
            -o $DIR/unitigs/\$name \
            $DIR/graphs/\$name.dbg > $DIR/graphs/\$name.clean.log 2>&1 \
        && rm $DIR/graphs/\$name.dbg*"


sbatch -J "build_joint_graph" \
     -o $DIR/logs/build_joint_graph.slog \
     -t 00-12 \
     --cpus-per-task 8 \
     --mem-per-cpu=10G \
     --partition=compute \
    --wrap="find $DIR/unitigs -name \"*.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
                -k 31 \
                --mode canonical \
                --mem-cap-gb 50 \
                --disk-swap ~/metagenome/scratch/nobackup \
                -p 8 \
                -o $DIR/graph; \
        /usr/bin/time -v $METAGRAPH transform -v \
                --to-fasta \
                --primary-kmers \
                -p 8 \
                -o $DIR/primary_contigs \
                $DIR/graph.dbg \
            && rm $DIR/graph.dbg; \
        /usr/bin/time -v $METAGRAPH build -v \
                -k 31 \
                --mode primary \
                --mem-cap-gb 50 \
                --disk-swap ~/metagenome/scratch/nobackup \
                -p 8 \
                -o $DIR/graph_primary \
                $DIR/primary_contigs.fasta.gz; \
        /usr/bin/time -v $METAGRAPH transform -v \
                --state small \
                -p 8 \
                -o $DIR/graph_primary_small \
                $DIR/graph_primary.dbg"



WINDOW_SIZE=5;
DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
mkdir ${DIR};
mkdir ${DIR}/logs;
mkdir ${DIR}/graphs;
mkdir ${DIR}/unitigs;

sbatch -J "build_clean_graphs" \
     --array=1-9858 \
     -o $DIR/logs/build_clean_graphs.slog \
     -t 00-12 \
     --cpus-per-task 4 \
     --mem-per-cpu=15G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="file=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${DIR}/../input_fastq.txt); \
        file=\${file%_1.fastq.gz}; \
        name=\$(basename \${file}); \
        [ -f $DIR/unitigs/\$name.fasta.gz ] && exit 0; \
        /usr/bin/time -v $METAGRAPH build -v -k 31 \
            --mode canonical \
            --count-kmers --count-width 16 \
            --mem-cap-gb 30 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 4 \
            -o $DIR/graphs/\$name \
            $DIR/../unitigs/\$name.fasta.gz > $DIR/graphs/\$name.build.log 2>&1 \
        && /usr/bin/time -v $METAGRAPH transform \
            --to-fasta --primary-kmers \
            --smoothing-window ${WINDOW_SIZE} \
            -p 4 \
            -o $DIR/unitigs/\$name \
            $DIR/graphs/\$name.dbg > $DIR/graphs/\$name.transform.log 2>&1 \
        && rm $DIR/graphs/\$name.dbg*"



DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
cd $DIR;
mkdir -p batches;
cd batches;
split -d -n r/40 <(find $DIR/unitigs -name "*.fasta.gz" | shuf);
mkdir -p ${DIR}/columns;

DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    sbatch -J "gtex_annotate_${WINDOW_SIZE}_${list}" \
         -o ${DIR}/logs/annotate_${list}.slog \
         -t 00-4 \
         --cpus-per-task 18 \
         --mem-per-cpu=2G \
        --wrap="cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/../graph_primary.dbg \
                --anno-filename \
                --separately \
                --count-kmers --count-width 16 \
                -o ${DIR}/columns \
                -p 18"; \
done


DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/../graph_primary.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
sbatch -J "gtex_count_${WINDOW_SIZE}_rd" \
     -o ${DIR}/logs/count_rd.lsf \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34 && \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34  && \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34";

DIR=~/metagenome/data/gtex_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
sbatch -J "gtex_count_${WINDOW_SIZE}_rd_brwt" \
     -d afterok:$(get_jobid gtex_count_${WINDOW_SIZE}_rd) \
     -o ${DIR}/logs/count_rd_brwt.slog \
     -t 1-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34 --parallel-nodes 10 && \
        /usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 34 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_int_brwt.annodbg";
```



```bash
################################## No counts ##################################

DIR=~/metagenome/data/gtex_counts/build/nobackup;
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/graph_primary.dbg ${DIR}/rd/graph.dbg;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

sbatch -J "gtex_rd_0" \
     -o $DIR/logs/rd_0.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34";

sbatch -J "gtex_rd_1" \
     -d afterok:$(get_jobid gtex_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34";

sbatch -J "gtex_rd_2" \
     -d afterok:$(get_jobid gtex_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34";

sbatch -J "gtex_rd_brwt" \
     -d afterok:$(get_jobid gtex_rd_2) \
     -o $DIR/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/rd/rd_columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34 --parallel-nodes 10";

sbatch -J "gtex_rd_brwt_relax" \
     -d afterok:$(get_jobid gtex_rd_brwt) \
     -o $DIR/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 34 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt.annodbg";
```
