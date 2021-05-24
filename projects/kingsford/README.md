```bash
find ~/metagenome/data/kingsford/kmc_20/ -name "*kmc_suf" > ~/metagenome/data/kingsford/kingsford_list_kmc.txt

bsub -J "build_single[1-2652]%500" \
     -o ~/metagenome/data/kingsford/build_single.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ~/metagenome/data/kingsford/kingsford_list_kmc.txt); \
        /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA build -v \
            -k 20 \
            --mode canonical \
            --count-kmers --count-width 32 \
            --mem-cap-gb 8 \
            -p 8 \
            -o ~/metagenome/data/kingsford/single_graphs/\\\$(basename \\\${file%.fasta.gz.kmc.kmc_suf}) \
            \\\$file \
        2>&1"

mkdir ~/metagenome/data/kingsford/smoothing_1000;
mkdir ~/metagenome/data/kingsford/smoothing_1000/logs;
mkdir ~/metagenome/data/kingsford/smoothing_1000/unitigs;

DIR=~/metagenome/data/kingsford/smoothing_1000;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "extract_contigs[1-2652]%500" \
     -o $DIR/logs/extract_contigs.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=2000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ~/metagenome/data/kingsford/kingsford_list_kmc.txt); \
        file=\\\$(basename \\\${file%.fasta.gz.kmc.kmc_suf})
        /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --primary-kmers \
            --smoothing-window 1000 \
            -p 8 \
            -o $DIR/unitigs/\\\$file \
            ~/metagenome/data/kingsford/single_graphs/\\\${file}.dbg \
        2>&1"

bsub -J "kingsford_build" \
     -oo $DIR/logs/build_canonical.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=3000] span[hosts=1]" \
    "find $DIR/unitigs -name \"*.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 20 \
            --mode canonical \
            --mem-cap-gb 80 \
            -p 72 \
            -o $DIR/kingsford_canonical \
            2>&1"

bsub -J "kingsford_to_primary" \
     -w "kingsford_build" \
     -oo $DIR/logs/extract_primary_contigs.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=2000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta --primary-kmers \
            -o $DIR/kingsford \
            $DIR/kingsford_canonical.dbg \
            -p 72 \
            2>&1"

bsub -J "kingsford_build" \
     -w "kingsford_to_primary" \
     -oo ~/metagenome/data/kingsford/build_primary.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1500] span[hosts=1]" \
    "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
            -k 20 \
            --mode primary \
            --mem-cap-gb 40 \
            -p 36 \
            -o ~/metagenome/data/kingsford/kingsford_primary \
            ~/metagenome/data/kingsford/kingsford_primary.fasta.gz \
            2>&1"


DIR=~/metagenome/data/kingsford/smoothing_1000;
cd $DIR
mkdir batches
cd batches
split -l 150 <(find $DIR/unitigs -name "*.fasta.gz" | shuf)
mkdir ${DIR}/columns

DIR=~/metagenome/data/kingsford/smoothing_1000;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd $DIR/batches;
for list in x*; do
    bsub -J "kingsford_annotate_${list}" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH annotate -v \
                -i ~/metagenome/data/kingsford/kingsford_primary.dbg \
                --anno-filename \
                --separately \
                -o ${DIR}/columns \
                -p 36 \
                2>&1 | tee ${DIR}/logs/annotate_${list}.log"; \
done


DIR=~/metagenome/data/kingsford/annotation;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd $DIR
bsub -J "kingsford_brwt" \
     -oo ${DIR}/logs/brwt.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1000] span[hosts=1]" \
    "find ${DIR}/columns_primary -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type brwt \
            --greedy --fast --subsample 1000000 \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/brwt.log"

bsub -J "kingsford_brwt_relax" \
     -w "kingsford_brwt" \
     -oo ${DIR}/logs/brwt_relax.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 72 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.brwt.annodbg \
        2>&1 | tee ${DIR}/logs/brwt_relax.log"
```





## Index with k-mer counts

```bash
WINDOW_SIZE=1;
mkdir ~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
mkdir ~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE}/logs;
mkdir ~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE}/unitigs;

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "extract_contigs[1-2652]%500" \
     -o $DIR/logs/extract_contigs.lsf \
     -W 4:00 \
     -n 4 -R "rusage[mem=2000] span[hosts=1]" \
        "file=\\\$(sed -n \${LSB_JOBINDEX}p ~/metagenome/data/kingsford/kingsford_list_kmc.txt); \
        file=\\\$(basename \\\${file%.fasta.gz.kmc.kmc_suf})
        /usr/bin/time -v $METAGRAPH clean \
            --to-fasta --primary-kmers \
            --smoothing-window ${WINDOW_SIZE} \
            -p 8 \
            -o $DIR/unitigs/\\\$file \
            ~/metagenome/data/kingsford/single_graphs/\\\${file}.dbg \
        2>&1"

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
bsub -J "split" \
     -w "extract_contigs" \
     -o /dev/null -W 1:00 -n 1 -R "rusage[mem=1000]" \
        "cd $DIR; \
        mkdir -p batches; \
        cd batches; \
        split -d -n r/40 <(find $DIR/unitigs -name "*.fasta.gz" | shuf); \
        mkdir -p ${DIR}/columns;";

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "kingsford_annotate_${list}" \
         -w "split" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i ~/metagenome/data/kingsford/kingsford_primary.dbg \
                --anno-filename \
                --separately \
                --count-kmers --count-width 32 \
                -o ${DIR}/columns \
                -p 36 \
                2>&1 | tee ${DIR}/logs/annotate_${list}.log"; \
done

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s ~/metagenome/data/kingsford/kingsford_primary.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
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

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
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

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
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

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "kingsford_count_rd_brwt" \
     -w "kingsford_count_rd_2" \
     -oo ${DIR}/logs/count_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/count_rd_brwt.log";

DIR=~/metagenome/data/kingsford/nobackup/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
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
