## Prepare and build graph

```bash
KMC=~/projects/projects2014-metagenome/metagraph/build_release/KMC/kmc;
DIR=~/metagenome/data/kingsford_21;
mkdir $DIR/kmc_21_filtered;
mkdir $DIR/logs;

for cutoff in {1,3,10,20,50}; do
  ids=$DIR/kingsford_${cutoff}.txt;
  bsub -J "filter[1-$(cat $ids | wc -l)]%100" \
       -o $DIR/logs/kmc_count.lsf \
       -W 4:00 \
       -n 1 -R "rusage[mem=8000] span[hosts=1]" \
            "id=\\\$(sed -n \${LSB_JOBINDEX}p $ids); \
            mkdir ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache; \
            file=~/metagenome/raw_data/kingsford/data_fasta/\\\${id}.fasta.gz; \
            /usr/bin/time -v $KMC -k21 -m6 -sm -ci$cutoff -fm -t2 \
                \\\${file} \
                $DIR/kmc_21_filtered/\\\$id \
                ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache \
            2>&1 | tee $DIR/kmc_21_filtered/\\\${id}.log;
            rm -r ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache"
done

find $DIR/kmc_21_filtered -name "*kmc_suf" > $DIR/kmc_list.txt

mkdir $DIR/unitigs;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_single[1-2652]%500" \
     -o $DIR/logs/build_single.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/kmc_list.txt); \
        /usr/bin/time -v $METAGRAPH build \
            -k 21 \
            --mode canonical \
            --mem-cap-gb 8 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 8 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            \\\$file \
        2>&1 | tee $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).log; \
        /usr/bin/time -v $METAGRAPH transform \
            --to-fasta --primary-kmers \
            -p 8 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg \
        2>&1 | tee -a $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).log; \
        rm $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg*"

bsub -J "build_graph" \
     -w "build_single" \
     -oo $DIR/logs/build_graph.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=3000] span[hosts=1]" \
    "find $DIR/unitigs -name \"*.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode canonical \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/kingsford_canonical \
            2>&1 | tee $DIR/logs/build_graph.log; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta --primary-kmers \
            -o $DIR/kingsford_primary \
            $DIR/kingsford_canonical.dbg \
            -p 72 \
            2>&1 | tee -a $DIR/logs/build_graph.log; \
    rm $DIR/kingsford_canonical.dbg; \
    /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode primary \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/kingsford \
            $DIR/kingsford_primary.fasta.gz \
            2>&1 | tee -a $DIR/logs/build_graph.log; \
    rm $DIR/kingsford_primary.fasta.gz"
```


## Index with k-mer counts

```bash
WINDOW_SIZE=1;
# WINDOW_SIZE=1000000000;
DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
mkdir $DIR;
mkdir $DIR/logs;
mkdir $DIR/unitigs;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_single_${WINDOW_SIZE}[1-2652]%500" \
     -o $DIR/logs/build_single.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/../kmc_list.txt); \
        /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode canonical \
            --count-kmers --count-width 32 \
            --mem-cap-gb 8 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 8 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            \\\$file \
        2>&1 | tee $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).log; \
        /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --primary-kmers \
            --smoothing-window ${WINDOW_SIZE} \
            -p 8 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg \
        2>&1 | tee -a $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).log; \
        rm $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg*"

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
bsub -J "split_${WINDOW_SIZE}" \
     -w "build_single_${WINDOW_SIZE}" \
     -o /dev/null -W 1:00 -n 1 -R "rusage[mem=1000]" \
        "cd $DIR; \
        mkdir -p batches; \
        cd batches; \
        split -d -n r/40 <(find $DIR/unitigs -name "*.fasta.gz" | shuf); \
        mkdir -p ${DIR}/columns;";

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "kingsford_annotate_${list}" \
         -w "split_${WINDOW_SIZE}" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/../kingsford.dbg \
                --anno-filename \
                --separately \
                --count-kmers --count-width 32 \
                -o ${DIR}/columns \
                -p 36 \
                2>&1 | tee ${DIR}/logs/annotate_${list}.log"; \
done

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_count_rd_0" \
     -w "kingsford_annotate_*" \
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

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_count_rd_1" \
     -w "kingsford_count_rd_0" \
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

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_count_rd_2" \
     -w "kingsford_count_rd_1" \
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

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_count_rd_brwt" \
     -w "kingsford_count_rd_2" \
     -oo ${DIR}/logs/count_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/count_rd_brwt.log";

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_count_rd_brwt_relax" \
     -w "kingsford_count_rd_brwt" \
     -oo ${DIR}/logs/count_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=5000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_int_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/count_rd_brwt_relax.log";
```

### Query

```bash
DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "kingsford_query_old" \
     -oo ${DIR}/logs/query_rd_brwt_relax_old.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=30000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH query --query-mode matches -v \
            -i ${DIR}/kingsford_small.dbg \
            -a ${DIR}/annotation_old.relaxed.row_diff_brwt.annodbg \
            ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa \
            2>&1 | tee ${DIR}/logs/query_rd_brwt_relax_old.log";
bsub -J "kingsford_query" \
     -oo ${DIR}/logs/query_rd_brwt_relax.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=30000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH query --query-mode matches -v \
            -i ${DIR}/kingsford_small.dbg \
            -a ${DIR}/annotation.relaxed.row_diff_brwt.annodbg \
            ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa \
            2>&1 | tee ${DIR}/logs/query_rd_brwt_relax.log";


WINDOW_SIZE=1;
# WINDOW_SIZE=1000000000;
DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

bsub -J "kingsford_count_query" \
     -oo ${DIR}/logs/query_count_rd_brwt_relax.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=30000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH query --query-mode counts-sum -v \
            -i ${DIR}/../kingsford_small.dbg \
            -a ${DIR}/annotation.relaxed.row_diff_int_brwt.annodbg \
            ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa \
            2>&1 | tee ${DIR}/logs/query_count_rd_brwt_relax.log";
```


## Binary annotation without k-mer counts

```bash
DIR=~/metagenome/data/kingsford_21;
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_0" \
     -w "kingsford_annotate_*" \
     -oo ${DIR}/logs/rd_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/rd_0.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_1" \
     -w "kingsford_rd_0" \
     -oo ${DIR}/logs/rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/rd_1.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_2" \
     -w "kingsford_rd_1" \
     -oo ${DIR}/logs/rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/rd_2.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_brwt" \
     -w "kingsford_rd_2" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/rd_brwt.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_brwt_relax" \
     -w "kingsford_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=5000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/rd_brwt_relax.log";
```


## RowDiff 1.0

```bash
git checkout 0d9feb76a9840b92031c25571cbc0f23ffd1cbe2

DIR=~/metagenome/data/kingsford_21;
mkdir $DIR/rd_old;
mkdir $DIR/rd_old/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd_old/graph.dbg;

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_1" \
     -w "kingsford_annotate_*" \
     -oo ${DIR}/logs/old_rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --max-path-length 100 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/rd_old/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/old_rd_1.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_2" \
     -w "kingsford_rd_old_1" \
     -oo ${DIR}/logs/old_rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --optimize \
            --max-path-length 100 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/rd_old/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/old_rd_2.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_brwt" \
     -w "kingsford_rd_old_2" \
     -oo ${DIR}/logs/old_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd_old/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/annotation_old \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/old_rd_brwt.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_brwt_relax" \
     -w "kingsford_rd_old_brwt" \
     -oo ${DIR}/logs/old_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=5000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation_old.relaxed \
            ${DIR}/annotation_old.row_diff_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/old_rd_brwt_relax.log";
```


## Index with k-mer coordinates (lossless read compression)

### Compress with tools for seq compression
```bash
DIR=~/metagenome/data/kingsford;
mkdir $DIR/compressed;
mkdir $DIR/compressed/logs;
find ~/metagenome/raw_data/kingsford/data_fasta -name "*.fasta.gz" > $DIR/compressed/list.txt;
list=$DIR/compressed/list.txt;
bsub -J "noheader_gzip[1-$(cat $list | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/noheader_gzip.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          /usr/bin/time -v zcat \\\$file | sed 's/^>.*/>/' | gzip -9 > $DIR/compressed/\\\${id}_no_header.fasta.gz"

bsub -J "spring[1-$(cat $list | wc -l)]" \
     -w "noheader_gzip[*]" \
     -o ~/metagenome/data/kingsford/compressed/logs/noheader_spring.lsf \
     -W 24:00 \
     -n 8 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          /usr/bin/time -v spring -c -g --fasta-input -t 16 \
                                  -i $DIR/compressed/\\\${id}_no_header.fasta.gz \
                                  -w ${DIR}/compressed \
                                  -o ${DIR}/compressed/\\\${id}_no_header.spring"


bsub -J "gzip[1-$(cat $list | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/gzip.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\$file); \
          /usr/bin/time -v zcat \\\$file | sed '/^>/d' | tr -d '\n' | gzip -9 > $DIR/compressed/\\\$id; \
          /usr/bin/time -v zcat \\\$file | sed '/^>/!d' | tr -d '\n' | gzip -9 > $DIR/compressed/\\\${id%.gz}_headers.gz"

bsub -J "spring[1-$(cat $list | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/spring.lsf \
     -W 24:00 \
     -n 8 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\$file); \
          /usr/bin/time -v spring -c -g --fasta-input -t 16 \
                                  -i \\\$file \
                                  -w ${DIR}/compressed \
                                  -o ${DIR}/compressed/\\\${id%.gz}.spring"
```

```bash
N=1
K=31
DIR=~/metagenome/data/kingsford_${K}_coordinates_dac_$N;
mkdir $DIR;
mkdir $DIR/logs;

cd $DIR
# cat ~/metagenome/data/kingsford_${K}/non_empty.txt | shuf | head -n $N > input_files.txt
cp ~/metagenome/data/kingsford_21_coordinates_$N/input_files.txt $DIR/

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "build_graph_${K}_${N}" \
     -oo $DIR/logs/build_graph.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat $DIR/input_files.txt \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k ${K} \
            --mem-cap-gb 100 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/kingsford \
            2>&1 | tee $DIR/logs/build_graph.log; \
    $METAGRAPH transform --state small -o $DIR/kingsford_small $DIR/kingsford.dbg -p 72"

mkdir -p ${DIR}/columns;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "annotate_${K}_${N}_[1-$(cat $DIR/input_files.txt | wc -l)]%100" \
     -w "build_graph_${K}_${N}" \
     -o ${DIR}/logs/annotate.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "file=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/input_files.txt); \
    id=\\\$(basename \\\$file); \
        /usr/bin/time -v $METAGRAPH annotate \
            -i $DIR/kingsford.dbg \
            --coordinates \
            --anno-filename \
            -o ${DIR}/columns/\\\$id \
            \\\$file \
            -p 72 \
            2>&1 | tee -a ${DIR}/logs/annotate.log"; \

mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/kingsford.dbg ${DIR}/rd/graph.dbg;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "kingsford_${K}_${N}_rd_0" \
     -w "annotate_${K}_${N}_*" \
     -oo ${DIR}/logs/coord_rd_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --coordinates \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/coord_rd_0.log";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "kingsford_${K}_${N}_rd_1" \
     -w "kingsford_${K}_${N}_rd_0" \
     -oo ${DIR}/logs/coord_rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --coordinates \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/coord_rd_1.log";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "kingsford_${K}_${N}_rd_2" \
     -w "kingsford_${K}_${N}_rd_1" \
     -oo ${DIR}/logs/coord_rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --coordinates \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/coord_rd_2.log";


K=31; N=1; echo Data: $(ls -l ../kingsford_${K}_coordinates_$N/compressed/*.gz | sizeb); echo Anno: $(( $(ls -l ../kingsford_${K}_coordinates_$N/rd/rd_columns/*.annodbg* | sizeb) + $(ls -l ../kingsford_${K}_coordinates_$N/rd/graph.dbg.anchors | sizeb) )); echo Graph: $(ls -l ../kingsford_${K}_coordinates_$N/kingsford_small.dbg | sizeb);

for N in {1,3,9}; do K=31; \
    DIR=../kingsford_${K}_coordinates_$N;
    echo N: $(ls $DIR/rd/rd_columns/*.coord* | wc -l); \
    echo Data: $(for x in $DIR/rd/rd_columns/*.coord*; do \
                    ls -l ../kingsford/compressed/$(basename ${x%.fasta.gz.column.annodbg.coords})_no_header.fasta.gz; \
                 done | sizeb); \
    echo Spring: $(for x in $DIR/rd/rd_columns/*.coord*; do \
                    ls -l ../kingsford/compressed/$(basename ${x%.fasta.gz.column.annodbg.coords})_no_header.spring; \
                 done | sizeb); \
    echo Anno: $(( $(ls -l $DIR/rd/rd_columns/*.annodbg* | sizeb) \
                    + $(ls -l $DIR/rd/graph.dbg.anchors | sizeb) )); \
    echo Graph: $(ls -l $DIR/kingsford_small.dbg | sizeb); echo ""; \
done

for N in {1,3,9,30,100}; do K=31; \
    arr=( )
    for x in ../kingsford_${K}_coordinates_$N/rd/rd_columns/*.coord*; do \
        L=$(zless ~/metagenome/raw_data/kingsford/data_fasta/$(basename ${x%.column.annodbg.coords}) | head -n 1 | grep -Eo '[0-9]+$'); \
        if [[ $L -ge 70 ]]; then
            arr+=($x)
        fi
    done
    echo N: $(for x in ${arr[@]}; do echo $x; done | wc -l); \
    echo Data: $(for x in ${arr[@]}; do \
                    ls -l ../kingsford/compressed/$(basename ${x%.fasta.gz.column.annodbg.coords})_no_header.spring; \
                 done | sizeb); \
    echo Anno: $(( $(for x in ${arr[@]}; do ls -l ${x%.coords}*; done | sizeb) \
                    + $(ls -l ../kingsford_${K}_coordinates_$N/rd/graph.dbg.anchors | sizeb) )); \
    echo Graph: $(ls -l ../kingsford_${K}_coordinates_$N/kingsford_small.dbg | sizeb); echo -e ""; \
done
```
