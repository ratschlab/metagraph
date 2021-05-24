```bash
mkdir ~/metagenome/data/gtex_counts/build/nobackup;
mkdir ~/metagenome/data/gtex_counts/build/nobackup/logs;
mkdir ~/metagenome/data/gtex_counts/build/nobackup/graphs;
mkdir ~/metagenome/data/gtex_counts/build/nobackup/unitigs;

find /cluster/work/grlab_shared/GTEx/extract/dbGaP-9608/fastq -name "*_1.fastq.gz" \
    > ~/metagenome/data/gtex_counts/build/nobackup/input_fastq.txt

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "build_graphs[1-9858]" \
     -o $DIR/logs/build_graphs.lsf \
     -W 12:00 \
     -n 8 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ${DIR}/input_fastq.txt); \
        file=\\\${file%_1.fastq.gz}; \
        name=\\\$(basename \\\${file}); \
        /usr/bin/time -v $METAGRAPH build -v -k 31 \
            --mode canonical \
            --count-kmers --count-width 12 \
            --mem-cap-gb 30 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 16 \
            -o $DIR/graphs/\\\$name \
            \\\${file}_1.fastq.gz \\\${file}_2.fastq.gz \
        2>&1"

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "extract_contigs[1-9858]" \
     -w "exit(build_graphs[*])" \
     -o $DIR/logs/extract_contigs.lsf \
     -W 12:00 \
     -n 8 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ${DIR}/input_fastq.txt); \
        file=\\\${file%_1.fastq.gz}; \
        name=\\\$(basename \\\${file}); \
        /usr/bin/time -v $METAGRAPH clean \
            --to-fasta --primary-kmers \
            --prune-tips 62 \
            --prune-unitigs 0 \
            --fallback 2 \
            -p 16 \
            -o $DIR/unitigs/\\\$name \
            $DIR/graphs/\\\$name.dbg \
        2>&1 && rm \"$DIR/graphs/\\\$name.*\""

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "build_graph" \
     -o $DIR/logs/build.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "find $DIR/unitigs -name \"*.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode canonical \
            --mem-cap-gb 100 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/graph \
        2>&1"

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "extract_primary" \
     -w "build_graph" \
     -o $DIR/logs/extract_primary_contigs.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta \
            --primary-kmers \
            -p 72 \
            -o $DIR/primary_contigs \
            $DIR/graph.dbg \
        2>&1"

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "build_primary" \
     -w "extract_primary" \
     -o $DIR/logs/build_primary.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode primary \
            --mem-cap-gb 100 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/graph_primary \
            $DIR/primary_contigs.fasta.gz \
        2>&1"

DIR=~/metagenome/data/gtex_counts/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
find $DIR/unitigs -name "*.fasta.gz" | shuf > ${DIR}/unitigs.txt;

bsub -J "build_clean_graphs[1-9759]" \
     -o $DIR/logs/build_clean_graphs.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ${DIR}/unitigs.txt); \
        file=\\\${file%.fasta.gz}; \
        name=\\\$(basename \\\${file}); \
        /usr/bin/time -v $METAGRAPH build -v -k 31 \
            --mode canonical \
            --count-kmers --count-width 12 \
            --mem-cap-gb 30 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 8 \
            -o $DIR/graphs/\\\$name \
            \\\${file}.fasta.gz \
        2>&1"

WINDOW_SIZE=31;
DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
mkdir ${DIR};
mkdir ${DIR}/logs;
mkdir ${DIR}/unitigs;

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "extract_contigs[1-9759]%500" \
     -o $DIR/logs/extract_contigs.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ${DIR}/../unitigs.txt); \
        file=\\\${file%.fasta.gz}; \
        name=\\\$(basename \\\${file}); \
        /usr/bin/time -v $METAGRAPH clean \
            --to-fasta --primary-kmers \
            --smoothing-window ${WINDOW_SIZE} \
            -p 8 \
            -o $DIR/unitigs/\\\$name \
            $DIR/../graphs/\\\${name}.dbg \
        2>&1"

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
cd $DIR;
mkdir -p batches;
cd batches;
split -d -n r/40 <(find $DIR/unitigs -name "*.fasta.gz" | shuf);
mkdir -p ${DIR}/columns;

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "gtex_annotate_${WINDOW_SIZE}_${list}" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/../graph_primary.dbg \
                --anno-filename \
                --separately \
                --count-kmers --count-width 12 \
                -o ${DIR}/columns \
                -p 36 \
                2>&1 | tee ${DIR}/logs/annotate_${list}.log"; \
done


DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/../graph_primary.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "gtex_count_${WINDOW_SIZE}_rd_0" \
     -w "gtex_annotate_${WINDOW_SIZE}_*" \
     -oo ${DIR}/logs/count_rd_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/count_rd_0.log";

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "gtex_count_${WINDOW_SIZE}_rd_1" \
     -w "gtex_count_${WINDOW_SIZE}_rd_0" \
     -oo ${DIR}/logs/count_rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/count_rd_1.log";

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "gtex_count_${WINDOW_SIZE}_rd_2" \
     -w "gtex_count_${WINDOW_SIZE}_rd_1" \
     -oo ${DIR}/logs/count_rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/count_rd_2.log";

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "gtex_count_${WINDOW_SIZE}_rd_brwt" \
     -w "gtex_count_${WINDOW_SIZE}_rd_2" \
     -oo ${DIR}/logs/count_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/count_rd_brwt.log";

DIR=~/metagenome/data/gtex_counts/build/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "gtex_count_${WINDOW_SIZE}_rd_brwt_relax" \
     -w "gtex_count_${WINDOW_SIZE}_rd_brwt" \
     -oo ${DIR}/logs/count_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 16 -R "rusage[mem=20000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_int_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/count_rd_brwt_relax.log";
```
