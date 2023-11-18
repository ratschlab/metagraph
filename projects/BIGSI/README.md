# All-microbial Index

[BIGSI](https://bigsi.readme.io/docs/)

## Get preprocessed data
Download clean mccortex graphs from [here](http://ftp.ebi.ac.uk/pub/software/bigsi/nat_biotech_2018/ctx/) and extract unitigs:
```bash
cat accessions.txt | xargs -n 1 -P 15 ./get_data.sh 2>&1 | tee log.txt
```

Accession IDs were extracted from the dumped [metadata](http://ftp.ebi.ac.uk/pub/software/bigsi/nat_biotech_2018/all-microbial-index/metadata) with this command:
```bash
export LC_ALL=C;
cat ~/Downloads/metadata | sed -e 's/SRR/\
SRR/g' | sed -e 's/ERR/\
ERR/g' | sed -e 's/DRR/\
DRR/g' | sed -n -e 's/^\([SED]RR[0-9]\{1,7\}\)r.*$/\1/p' | sort | uniq > accessions.txt
```

Partition the files based on the IDs and move them to final destination
```bash
for ID in $(cat accessions.txt); do echo ${ID:0:6}; done | sort | uniq > list.txt

for GROUP in $(cat list.txt); do \
  mv $(find ~/metagenome/data/BIGSI/dumps/ -name "${GROUP}*.lsf") ~/metagenome/data/BIGSI/$GROUP/; \
done
```

## Build graph
```bash
for F in {\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}; do \
    bsub -J "build_${F}" \
         -oo ~/metagenome/data/BIGSI/graph.${F}.lsf \
         -W 24:00 \
         -n 15 -R "rusage[mem=23000] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/ -name \"*fasta.gz\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph build -v \
                -k 31 \
                --mode canonical \
                --parallel 30 \
                --mem-cap-gb 300 \
                --suffix $F \
                -o ~/metagenome/data/BIGSI/graph \
                2>&1 | tee ~/metagenome/data/BIGSI/graph.$F.log"; \
done

bsub -J "transform_${F}" -o /dev/null -W 8:00 -n 1 -R "rusage[mem=260000] span[hosts=1]" \
    "~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform \
        --state small -o ~/metagenome/data/BIGSI/graph ~/metagenome/data/BIGSI/graph_stat.dbg"
```

## Annotate
```bash
mkdir temp
cd temp
find ~/metagenome/data/BIGSI/data/ -name "*fasta.gz" > files_to_annotate.txt
split -l 5000 files_to_annotate.txt
cd ..

find ~/metagenome/data/BIGSI/data/ -name "*.column.annodbg"

cd temp
for list in x*; do
    bsub -J "annotate_${list}" \
         -oo ~/metagenome/data/BIGSI/logs/annotate_${list}.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=9000] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                -i ~/metagenome/data/BIGSI/graph.dbg \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/BIGSI/annotation/columns \
                2>&1"; \
done
cd ..

bsub -J "cluster" \
     -oo ~/metagenome/data/BIGSI/logs/cluster_columns.lsf \
     -W 120:00 \
     -n 48 -R "rusage[mem=42500] span[hosts=1]" \
    "find ~/metagenome/data/BIGSI/annotation/columns_canonical/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph transform_anno -v \
            --linkage --greedy \
            --subsample 5000000 \
            -o ~/metagenome/data/BIGSI/annotation/linkage_BIGSI.csv \
            --parallel 96 \
            2>&1";

bsub -J "cluster" \
     -oo ~/metagenome/data/BIGSI/logs/cluster_columns_1M.lsf \
     -W 120:00 \
     -n 48 -R "rusage[mem=37500] span[hosts=1]" \
    "find ~/metagenome/data/BIGSI/annotation/columns_canonical/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph transform_anno -v \
            --linkage --greedy \
            --subsample 1000000 \
            -o ~/metagenome/data/BIGSI/annotation/linkage_BIGSI_1M.csv \
            --parallel 96 \
            2>&1";

bsub -J "build_brwt_5M" \
     -oo ~/metagenome/data/BIGSI/logs/build_brwt_5M.lsf \
     -W 120:00 \
     -n 48 -R "rusage[mem=42000] span[hosts=1]" \
    "cat ~/metagenome/data/BIGSI/annotation/columns.txt \
        | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --anno-type brwt \
            --linkage-file ~/metagenome/data/BIGSI/annotation/linkage_BIGSI.csv \
            -o ~/metagenome/data/BIGSI/annotation/annotation_5M \
            -p 96 --parallel-nodes 48 \
            2>&1"

bsub -J "relax_brwt_5M" \
     -oo ~/metagenome/data/BIGSI/logs/relax_brwt_5M.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=17000] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
        -p 36 \
        --relax-arity 16 \
        -o /cluster/home/mikhaika/metagenome/data/BIGSI/annotation/annotation_5M.relaxed \
        /cluster/home/mikhaika/metagenome/data/BIGSI/annotation/annotation_5M.brwt.annodbg \
        2>&1"
```

## Transform to primary
```bash
bsub -J "BIGSI_to_primary" \
     -oo ~/metagenome/data/BIGSI/logs/transform_to_primary.lsf \
     -W 12:00 \
     -n 18 -R "rusage[mem=5600] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform -v \
            --to-fasta --primary-kmers \
            -o ~/metagenome/data/BIGSI/graph_primary \
            ~/metagenome/data/BIGSI/graph.dbg \
            -p 18 \
            2>&1"

bsub -J "BIGSI_primary_build" \
     -oo ~/metagenome/data/BIGSI/logs/build_primary.lsf \
     -W 12:00 \
     -n 36 -R "rusage[mem=6000] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
            -k 31 \
            --mode primary \
            --mem-cap-gb 100 \
            -o ~/metagenome/data/BIGSI/graph_primary \
            ~/metagenome/data/BIGSI/graph_primary.fasta.gz \
            --disk-swap ~/metagenome/scratch/ \
            --disk-cap-gb 1000 \
            -p 36 \
            2>&1"


mkdir temp
cd temp
find ~/metagenome/data/BIGSI/data/ -name "*fasta.gz" > files_to_annotate.txt
split -l 5000 files_to_annotate.txt
cd ..

cd temp
for list in x*; do
    bsub -J "annotate_${list}_primary" \
         -oo ~/metagenome/data/BIGSI/logs/annotate_${list}_primary.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=5000] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                -i ~/metagenome/data/BIGSI/graph_primary.dbg \
                --fwd-and-reverse \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/BIGSI/annotation/columns_primary \
                2>&1"; \
done
cd ..


find ~/metagenome/data/BIGSI/annotation/columns_primary/ -name "*.column.annodbg" > columns_primary.txt

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/BIGSI/build;
mkdir $DIR/logs;

sbatch -J "BIGSI_build" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/../data -name \"*.unitigs.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --inplace \
            --mode canonical \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup \
            -p 34 \
            -o $DIR/graph; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta \
            --primary-kmers \
            -p 34 \
            -o $DIR/primary_contigs \
            $DIR/graph.dbg && \
    rm $DIR/graph.dbg; \
    /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode primary \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup \
            -p 34 \
            -o $DIR/graph_primary \
            $DIR/primary_contigs.fasta.gz; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 34 \
            -o $DIR/graph_primary_small \
            $DIR/graph_primary.dbg;"


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/BIGSI/annotation;
mkdir ${DIR}/rd;
ln -s $DIR/../graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "BIGSI_rd_0" \
     -o $DIR/../logs/rd_0.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=20G \
     --exclude compute-biomed-10 \
    --wrap="cat $DIR/columns_primary.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 500 \
            -p 34"

sbatch -J "BIGSI_rd_1" \
     -d afterok:$(get_jobid BIGSI_rd_0) \
     -o $DIR/../logs/rd_1.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=20G \
     --exclude compute-biomed-10 \
    --wrap="cat $DIR/columns_primary.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 500 \
            -p 34 && \
        rm ${DIR}/rd/rd_columns/*.row_count"

sbatch -J "BIGSI_rd_2" \
     -d afterok:$(get_jobid BIGSI_rd_1) \
     -o $DIR/../logs/rd_2.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=20G \
     --exclude compute-biomed-10 \
    --wrap="cat $DIR/columns_primary.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 500 \
            -p 34 && \
        rm ${DIR}/rd/rd_columns/*.row_reduction && \
        rm ${DIR}/rd/graph.dbg.pred* && \
        rm ${DIR}/rd/graph.dbg.succ*"


sed 's/annotation\/columns_primary/annotation\/rd\/rd_columns/' $DIR/columns_primary.txt \
    | sed 's/.column.annodbg/.row_diff.annodbg/' > $DIR/rd/columns_rd.txt


sbatch -J "BIGSI_rd_brwt" \
     -d afterok:$(get_jobid BIGSI_rd_2) \
     -o $DIR/../logs/rd_brwt.slog \
     -t 5-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=18G \
    --wrap="cat $DIR/rd/columns_rd.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --linkage-file $DIR/linkage_BIGSI_primary.csv \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34 --parallel-nodes 10"

sbatch -J "BIGSI_rd_brwt_relax" \
     -d afterok:$(get_jobid BIGSI_rd_brwt) \
     -o $DIR/../logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=18G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"


# bsub -J "build_rbbrwt_relax20_primary" \
#      -oo ~/metagenome/data/BIGSI/logs/build_rbbrwt_primary.lsf \
#      -W 120:00 \
#      -n 36 -R "rusage[mem=13000] span[hosts=1]" \
#     "cat ~/metagenome/data/BIGSI/annotation/columns_primary.txt \
#         | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
#             --anno-type rb_brwt --relax-arity 20 \
#             -o ~/metagenome/data/BIGSI/annotation/annotation_primary \
#             -p 36 \
#             2>&1"


bsub -J "cluster" \
     -oo ~/metagenome/data/BIGSI/logs/cluster_columns_primary.lsf \
     -W 360:00 \
     -n 48 -R "rusage[mem=42500] span[hosts=1]" \
    "cat ~/metagenome/data/BIGSI/annotation/columns_primary.txt \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --linkage --greedy \
            --subsample 5000000 \
            -o ~/metagenome/data/BIGSI/annotation/linkage_BIGSI_primary.csv \
            --parallel 96 \
            2>&1";

bsub -J "build_brwt_primary" \
     -oo ~/metagenome/data/BIGSI/logs/build_brwt_primary.lsf \
     -W 72:00 \
     -n 48 -R "rusage[mem=20000] span[hosts=1]" \
    "cat ~/metagenome/data/BIGSI/annotation/columns_primary.txt \
        | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --anno-type brwt \
            --linkage-file ~/metagenome/data/BIGSI/annotation/linkage_BIGSI_primary.csv \
            -o ~/metagenome/data/BIGSI/annotation/annotation_primary \
            -p 48 --parallel-nodes 48 \
            2>&1"

bsub -J "relax_brwt_primary" \
     -oo ~/metagenome/data/BIGSI/logs/relax_brwt_primary.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=17000] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
        -p 36 \
        --relax-arity 16 \
        -o /cluster/home/mikhaika/metagenome/data/BIGSI/annotation/annotation_primary.relaxed \
        /cluster/home/mikhaika/metagenome/data/BIGSI/annotation/annotation_primary.brwt.annodbg \
        2>&1"
```

## Generate subsets
```bash
cat temp/files_to_annotate.txt | shuf > subsets/files_shuffled.txt
cd subsets
for i in {1..33}; do N=$((750 * i)); head -n $N files_shuffled.txt > files_$N.txt; done
```

### Build graphs for subsets
```bash
mkdir ~/metagenome/data/BIGSI/subsets
for i in {1..33}; do
    N=$((750 * i));
    bsub -J "subset_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.lsf \
         -W 24:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "cat subsets/files_${N}.txt \
            | /usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA build -v \
                -k 31 \
                --mode canonical \
                --parallel 30 \
                --mem-cap-gb 300 \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N} \
                2>&1"; \
done
```

### Transform graphs
```bash
for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_primary_contigs_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/to_primary_contigs_${N}.lsf \
         -W 12:00 \
         -n 20 -R "rusage[mem=2500] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform -v \
                --to-fasta --primary-kmers -p 20 \
                -o ~/metagenome/data/BIGSI/subsets/contigs/graph_subset_${N}.primary_contigs \
                ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.dbg \
                2>&1"; \
done


for i in {1..33}; do
    N=$((750 * i));
    bsub -J "subset_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/graph_subset_${N}_primary.lsf \
         -W 24:00 \
         -n 15 -R "rusage[mem=12000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA build -v \
                -k 31 \
                --mode primary \
                --parallel 30 \
                --mem-cap-gb 150 \
                --index-ranges 12 \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary \
                ~/metagenome/data/BIGSI/subsets/contigs/graph_subset_${N}.primary_contigs.fasta.gz \
                2>&1"; \
done


for i in {1..33}; do
    N=$((750 * i));
    bsub -J "to_small_graph_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/to_small_graph_${N}.lsf \
         -W 1:00 \
         -n 1 -R "rusage[mem=50000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph_DNA transform -v \
                --state small \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary.small.dbg \
                ~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary.dbg \
                2>&1"; \
done


for i in {1..33}; do
    N=$((750 * i));
    OUTDIR=~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary
    mkdir -p $OUTDIR
    bsub -J "annotate_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/annotate_subset_${N}_primary.lsf \
         -W 12:00 \
         -n 15 -R "rusage[mem=3000] span[hosts=1]" \
        "cat subsets/files_${N}.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA annotate -v \
                -i ~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary.dbg \
                --fwd-and-reverse \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o $OUTDIR \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_brwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_brwt_${N}_primary.lsf \
         -W 24:00 \
         -n 20 -R "rusage[mem=$((N + 1))] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_master/metagraph_DNA transform_anno -v \
                --anno-type brwt --greedy \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                --parallel 20 \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "relax_brwt_${N}" \
         -w "to_brwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/relax_brwt_${N}_primary.lsf \
         -W 24:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_master/metagraph_DNA relax_brwt -v \
                --relax-arity 32 \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary.relaxed \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary.brwt.annodbg \
                --parallel 20 \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_rbbrwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rbbrwt_${N}_primary.lsf \
         -W 24:00 \
         -n 20 -R "rusage[mem=$((7 * N / 4 + 1000))] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_master/metagraph_DNA transform_anno -v \
                --anno-type rb_brwt --relax-arity 32 \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                --parallel 20 \
                2>&1"; \
done
```


```bash
for i in {33..1}; do
    N=$((750 * i));
    mkdir ~/metagenome/data/BIGSI/subsets/graph_subset_${N}
    bsub -J "annotate_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/annotate_subset_${N}.lsf \
         -W 12:00 \
         -n 15 -R "rusage[mem=3000] span[hosts=1]" \
        "cat subsets/files_${N}.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA annotate -v \
                -i ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.dbg \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N} \
                2>&1"; \
done

for i in {33..1}; do
    N=$((750 * i));
    bsub -J "merge_columns_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/merge_columns_subset_${N}.lsf \
         -W 4:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/subsets/graph_subset_${N}/ -name \"*.annodbg\" | sort \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA merge_anno -v \
                -o ~/metagenome/data/BIGSI/subsets/annotation_subset_${N} \
                --parallel 20 \
                2>&1"; \
done
```

### Rename columns
```bash
./metagraph stats --print-col-names -a ~/metagenome/data/BIGSI/subsets/annotation_subset_24750.column.annodbg > ~/metagenome/data/BIGSI/subsets/rename_columns.txt
tail -n +3 ~/metagenome/data/BIGSI/subsets/rename_columns.txt > ~/metagenome/data/BIGSI/subsets/rename_columns_.txt
for x in $(cat ~/metagenome/data/BIGSI/subsets/rename_columns_.txt); do echo "$x $(basename ${x%.unitigs.fasta.gz})"; done > ~/metagenome/data/BIGSI/subsets/rename_columns.txt
rm ~/metagenome/data/BIGSI/subsets/rename_columns_.txt

extension="_primary.row_diff_sparse.annodbg";
for i in {33..1}; do
    N=$((750 * i));
    annotation=~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}${extension};
    new="${annotation%$extension}.new${extension}";
    bsub -J "rename_columns_${N}" \
         -oo /dev/null \
         -W 1:00 \
         -n 1 -R "rusage[mem=$((N * 6))] span[hosts=1]" \
        "~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno \
            --rename-cols ~/metagenome/data/BIGSI/subsets/rename_columns.txt \
            -o $new \
            $annotation \
            && mv $new $annotation"; \
done
```

## Transform annotation
```bash
for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_row_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_row_${N}_primary.lsf \
         -W 24:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph_DNA transform_anno -v \
                --anno-type row \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                --parallel 20 \
                2>&1"; \
done


for i in {1,5,10,33}; do
    N=$((750 * i));
    bsub -J "to_flat_${N}" \
         -w "to_row_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/row_to_flat_${N}_primary.lsf \
         -W 24:00 \
         -n 1 -R "rusage[mem=$((N * 10 + 15000))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type flat \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary.row.annodbg \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_rbfish_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/row_to_rbfish_${N}_primary.lsf \
         -W 20:00 \
         -n 1 -R "rusage[mem=$((N * 17))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type rbfish \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary.row.annodbg \
                2>&1"; \
done


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph_DNA;
for i in {1,5,10,33}; do
    N=$((750 * i));
    bsub -J "to_row_sparse_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/row_to_row_sparse_${N}_primary.lsf \
         -W 20:00 \
         -n 1 -R "rusage[mem=$((N * 17))] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_sparse \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary.row.annodbg \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_brwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_brwt_${N}.lsf \
         -W 24:00 \
         -n 20 -R "rusage[mem=$((N + 1))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type brwt --greedy \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.column.annodbg \
                --parallel 20 \
                2>&1"; \
done

for i in {33..1}; do
    N=$((750 * i));
    bsub -J "relax_brwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/relax_brwt_${N}.lsf \
         -W 96:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA relax_brwt -v \
                --relax-arity 20 \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.relaxed \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.brwt.annodbg \
                --parallel 20 \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_rbbrwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rbbrwt_${N}.lsf \
         -W 96:00 \
         -n 10 -R "rusage[mem=$((14 * N / 9 + 1000))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type rb_brwt --greedy \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.column.annodbg \
                --parallel 20 \
                2>&1"; \
done


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph_DNA;
for i in {33..1}; do
    N=$((750 * i));
    DIR=~/metagenome/data/BIGSI/subsets/annotation/rd_${N};
    mkdir $DIR;
    mkdir $DIR/rd_columns;
    ln -s ~/metagenome/data/BIGSI/subsets/graph_subset_${N}_primary.dbg $DIR/graph.dbg;
    sbatch -J "to_rd_0_${N}" \
         -o ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rd_0_${N}_primary.slog \
         -t 00-24 \
         --cpus-per-task 10 \
         --mem-per-cpu=12G \
        --wrap="find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.column.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 0 \
                --mem-cap-gb 100 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i $DIR/graph.dbg \
                -o $DIR/rd_columns/out \
                -p 20"; \
    sbatch -J "to_rd_1_${N}" \
         -d afterok:$(get_jobid to_rd_0_${N}) \
         -o ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rd_1_${N}_primary.slog \
         -t 00-24 \
         --cpus-per-task 10 \
         --mem-per-cpu=12G \
        --wrap="find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.column.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 1 \
                --mem-cap-gb 100 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i $DIR/graph.dbg \
                -o $DIR/rd_columns/out \
                -p 20"; \
    sbatch -J "to_rd_2_${N}" \
         -d afterok:$(get_jobid to_rd_1_${N}) \
         -o ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rd_2_${N}_primary.slog \
         -t 00-24 \
         --cpus-per-task 10 \
         --mem-per-cpu=12G \
        --wrap="find ~/metagenome/data/BIGSI/subsets/annotation/columns/graph_subset_${N}_primary/ -name \"*.column.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 2 \
                --mem-cap-gb 100 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i $DIR/graph.dbg \
                -o $DIR/rd_columns/out \
                -p 20"; \
    sbatch -J "to_rd_flat_${N}" \
         -d afterok:$(get_jobid to_rd_2_${N}) \
         -o ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rd_flat_${N}_primary.slog \
         -t 00-24 \
         --cpus-per-task 10 \
         --mem-per-cpu=12G \
        --wrap="find ${DIR}/rd_columns -name \"*.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_flat \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}_primary \
                -p 20"; \
done
```

## Query graph
```bash
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph

# file to query
for QUERY in ~/metagenome/data/BIGSI/subsets/query/samples/haib18CEM5453_HMCMJCCXY_SL336225.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/nucleotide_fasta_protein_homolog_model.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/SRR322873.fastq.gz \
                ~/metagenome/data/BIGSI/subsets/query/samples/DRR067889.fasta; do
    NAME=metagraph.rd_brwt.align_labeled
    # name of the output folder
    OUTDIR=~/metagenome/data/BIGSI/subsets/query_results/$(basename $QUERY)/${NAME}
    rm -rf $OUTDIR
    mkdir -p $OUTDIR

    for num_columns in $(seq 750 3000 24750); do
        # run="$METAGRAPH query -v --min-kmers-fraction-label 0.0 --query-mode matches \
        run="$METAGRAPH align -v \
                -i \${TMPDIR}/graph.dbg \
                -a \${TMPDIR}/graph.row_diff_brwt.annodbg \
                $QUERY"

        bsub -J "${NAME}.${num_columns}" \
            -W 24:00 \
            -n 36 -R "rusage[mem=4000] span[hosts=1] select[model==XeonGold_6140]" \
            -oo ${OUTDIR}/${num_columns}.lsf \
            " \
                TMPDIR=/dev/shm/$NAME_$num_columns_\${LSB_JOBID}; \
                mkdir \${TMPDIR}; \
                cp ~/metagenome/data/BIGSI/subsets/graph_subset_${num_columns}_primary.dbg \
                    \${TMPDIR}/graph.dbg; \
                cp ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${num_columns}_primary.relaxed.row_diff_brwt.annodbg \
                    \${TMPDIR}/graph.row_diff_brwt.annodbg; \
                /usr/bin/time -v $run > /dev/null 2> /dev/null; \
                for i in {1..3}; do \
                    /usr/bin/time -v $run > \${TMPDIR}/out 2>> \${TMPDIR}/err;
                done; \
                mv \${TMPDIR}/out ${OUTDIR}/${num_columns}.out; \
                mv \${TMPDIR}/err ${OUTDIR}/${num_columns}.err; \
                rm -r \${TMPDIR}"
    done
done
```

## Generate BIGSI bloom filters
```bash
bsub -J "bigsi_bloom[1-24750]%300" \
     -o ~/metagenome/data/BIGSI/subsets/bigsi/build_bloom.lsf \
     -W 8:00 \
     -n 1 -R "rusage[mem=3000] span[hosts=1]" \
     "file=\"\$(sed -n \${LSB_JOBINDEX}p subsets/files_24750.txt)\"; \
        x=\"\$(echo \$file | xargs -n 1 basename)\"; \
        sra_id=\"\${x%.unitigs.fasta.gz}\"; \
        ctx_file=\"/cluster/work/grlab/projects/metagenome/raw_data/BIGSI/ctx_subsets/\${sra_id}.ctx\"; \
        /usr/bin/time -v bigsi bloom \
            -c ~/projects/projects2014-metagenome/projects/BIGSI/bigsi_config.yaml \
            \$ctx_file \
            /cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/bigsi/bloom/\${sra_id}.bloom \
            2>&1"
```

```bash
mkdir bigsi_subsets

for i in {1..33}; do
    N=$((750 * i));
    file="bigsi_subsets/files_${N}.txt";
    rm -f $file;
    for f in $(cat subsets/files_${N}.txt); do
        sra_id="$(basename $f)";
        sra_id="${sra_id%.unitigs.fasta.gz}";
        echo -e "/cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/bigsi/bloom/${sra_id}.bloom\t$sra_id" >> $file;
    done
done


for i in {1..33}; do
    N=$((750 * i));
    mkdir ~/metagenome/data/BIGSI/subsets/bigsi/subsets/index_${N};
    cd ~/metagenome/data/BIGSI/subsets/bigsi/subsets/index_${N};
    bsub -J "bigsi_to_db_${N}" \
         -oo merge_blooms_${N}.lsf \
         -W $((N / 100)):00 \
         -n 5 -R "rusage[mem=$((N * 11))] span[hosts=1]" \
        "/usr/bin/time -v ~/anaconda3/bin/bigsi build \
                -f ~/projects/projects2014-metagenome/projects/BIGSI/bigsi_subsets/files_${N}.txt \
                -c ~/projects/projects2014-metagenome/projects/BIGSI/bigsi_config.yaml \
                2>&1"; \
done
```

## Generate COBS index
```bash
for i in {1..33}; do
    N=$((750 * i));
    mkdir -p "/cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/cobs/data/subset_${N}"
    for f in $(cat subsets/files_${N}.txt); do
        filename="$(basename $f)";
        ln -s "$f" "/cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/cobs/data/subset_${N}/${filename}";
    done
done

num_hashes=3;
fpr=10;
# num_hashes=4;
# fpr=5;
# num_hashes=7;
# fpr=1;
for i in {1..33}; do
    N=$((750 * i));
    bsub -J "cobs_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/cobs/build_cobs_compact_h${num_hashes}_fpr${fpr}_${N}.lsf \
         -W 24:00 \
         -n 10 -R "rusage[mem=11000] span[hosts=1]" \
        "/usr/bin/time -v ~/stuff/cobs/build/src/cobs compact-construct \
                -k 31 -c -f $(echo ${fpr}*0.01 | bc) -T 10 --num-hashes ${num_hashes} \
                /cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/cobs/data/subset_${N}/ \
                /cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/cobs/cobs_index_h${num_hashes}_fpr${fpr}_${N}.cobs_compact \
                2>&1"; \
done
```

#### Profile
```bash
#QUERY=~/metagenome/data/BIGSI/subsets/query/samples/haib18CEM5453_HMCMJCCXY_SL336225.fasta
QUERY=~/metagenome/data/BIGSI/subsets/query/samples/nucleotide_fasta_protein_homolog_model.fasta
#QUERY=~/metagenome/data/BIGSI/subsets/query/samples/DRR067889.fasta
./metagraph query -v --min-kmers-fraction-label 0.0 --query-mode matches \
    -i ~/metagenome/data/BIGSI/subsets/graph_subset_24750.dbg \
    -a ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_24750.rb_brwt.annodbg \
    $QUERY \
    > ~/metagenome/scratch/out_test.txt \
    2>&1
```