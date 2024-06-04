```bash
rm -rf ~/metagenome/data/refseq/release97/output_k31_annotation/rd;
mkdir ~/metagenome/data/refseq/release97/output_k31_annotation/rd;
mkdir ~/metagenome/data/refseq/release97/output_k31_annotation/rd/rd_columns;
ln -s ~/metagenome/data/refseq/release97/output_k31_graph_chunked/graph_k31.dbg \
      ~/metagenome/data/refseq/release97/output_k31_annotation/build/rd/graph_k31.dbg;
ln -s ~/metagenome/data/refseq/release97/output_k31_graph_chunked/graph_k31.edgemask \
      ~/metagenome/data/refseq/release97/output_k31_annotation/build/rd/graph_k31.edgemask;

DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
pushd ${DIR}/build/batches/input
ls -1S ~/metagenome/raw_data/refseq/release97/fna/complete \
    > ${DIR}/build/batches/input/input.txt
sed -i 's/^/\/cluster\/work\/grlab\/projects\/metagenome\/raw_data\/refseq\/release97\/fna\/complete\//' \
    ${DIR}/build/batches/input/input.txt
split input.txt -l 200
cat xae xai xam >> xaa
cat xaf xaj xan >> xab
cat xag xak xao xaq >> xac
cat xah xal xap xar >> xad
rm xah xal xap xar xag xak xao xaq xaf xaj xan xae xai xam
popd

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
cd ${DIR}/build/batches/input;
for list in x*; do
    bsub -J "annotate_RefSeq_${list}" \
         -oo ${DIR}/build/logs/annotate_${list}.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=20000] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH annotate -v \
                -i ${DIR}/graph.dbg \
                --parallel 72 \
                --cache 1 \
                --anno-header \
                --separately \
                -o ${DIR}/build/columns \
                2>&1 | tee ${DIR}/build/logs/annotate_${list}.log"; \
done

DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
mkdir ${DIR}/build/rd;
mkdir ${DIR}/build/rd/rd_columns;
ln -s ${DIR}/graph.dbg \
      ${DIR}/build/rd/graph.dbg;

cp ${DIR}/build/batches/input/input.txt ${DIR}/build/columns.txt
sed -i 's/raw_data\/refseq\/release97\/fna\/complete/data\/refseq\/release97\/output_k31_annotation\/k31_full\/build\/columns/' \
            ${DIR}/build/columns.txt
sed -i 's/$/.column.annodbg/' ${DIR}/build/columns.txt

DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
cd ${DIR}/build/batches/;
cat ../columns.txt | shuf > columns.txt;
split columns.txt -l 400


DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd ${DIR}/build/batches;
for list in x*; do
    bsub -J "refseq_rd_${list}" \
         -oo ${DIR}/build/logs/column_to_rd_${list}.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=20000] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${DIR}/build/batches/${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length 100 \
                --mem-cap-gb 400 \
                --parallel 72 \
                -o ${DIR}/build/rd/rd_columns/out \
                -i ${DIR}/build/rd/graph.dbg \
                2>&1 | tee ${DIR}/build/logs/column_to_rd_${list}.log"; \
done



DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd ${DIR}/build/batches;
for list in x*; do
    bsub -J "refseq_rd_${list}_opt" \
         -oo ${DIR}/build/logs/column_to_rd_${list}_opt.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${DIR}/build/batches/${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length 100 \
                --mem-cap-gb 400 \
                --parallel 72 \
                --optimize \
                -o ${DIR}/build/rd/rd_columns_opt/out \
                -i ${DIR}/build/rd/graph.dbg \
                2>&1 | tee ${DIR}/build/logs/column_to_rd_${list}_opt.log"; \
done


DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full/nobackup;
pushd ${DIR}/build/rd
cat ../col_names.txt | grep annotation | cut -d' ' -f6 > rd_columns.txt
sed -i "s/'//" rd_columns.txt
sed -i "s/'//" rd_columns.txt
sed -i 's/fna\/complete_split/output_k31_annotation\/k31_full\/nobackup\/build\/columns/' \
    columns.txt
popd


DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd ${DIR}/build/rd;
bsub -J "refseq_rd_brwt" \
     -oo ${DIR}/build/logs/build_rd_brwt.lsf \
     -W 72:00 \
     -n 48 -R "rusage[mem=43000] span[hosts=1]" \
    "cat ${DIR}/build/rd/rd_columns.txt \
        | gtime -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            -i ${DIR}/build/rd/graph.dbg \
            --linkage-file ${DIR}/build/columns.taxid.linkage.txt \
            -o ${DIR}/annotation \
            --disk-swap \\\"\\\" \
            -p 72 --parallel-nodes 15 \
            2>&1 | tee ${DIR}/build/logs/build_rd_brwt.log"


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_full/nobackup;
bsub -J "refseq_rd_brwt_relax_32" \
     -oo ${DIR}/build/logs/build_rd_brwt_relax.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=25000] span[hosts=1]" \
    "gtime -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed_32 \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/build/logs/build_rd_brwt_relax_32.log"


### BY TAXID



DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
mkdir -p ${DIR}/batches/input
pushd ${DIR}/batches/input
ls -1S ~/metagenome/raw_data/refseq/release97/fna/complete_split \
    > ${DIR}/batches/input/input.txt
sed -i 's/^/\/cluster\/work\/grlab\/projects\/metagenome\/raw_data\/refseq\/release97\/fna\/complete_split\//' \
    ${DIR}/batches/input/input.txt
for i in {1..6}; do
    sed -n "${i}~6p" ${DIR}/batches/input/input.txt > ${DIR}/batches/input/batch_${i}.txt;
done
popd

ln -s /cluster/work/grlab/projects/metagenome/data/refseq/release97/output_k31_graph_chunked/graph_k31.dbg \
        ${DIR}/rd/graph.dbg

cat ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/batches/input/input.txt \
        | xargs -P 200 -I {} sh -c "zcat {} | gzip -9 > {}.9"
cat ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/batches/input/input.txt \
        | xargs -P 200 -I {} sh -c "zcat {} | grep '>' >> ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/headers.txt"
cat ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/batches/input/input.txt | xargs -P 200 -I {} sh -c "zcat {} | grep -v '>' | tr -d '\n' | wc -c" | awk "{sum+=\$1}END{print sum}"
cat ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/batches/input/input.txt \
        | xargs -P 200 -I {} sh -c "ls -l {}.9" | sizeb
cat ~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build/logs/rd_2*.lsf | grep "reduced from" | cut -d' ' -f12 | awk "{sum+=\$1}END{print sum}"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
cd ${DIR}/batches/input;
#for list in {batch_1.txt,batch_2.txt,batch_3.txt,batch_4.txt,batch_5.txt,batch_6.txt}; do
for list in {batch_3.txt,batch_4.txt,batch_5.txt,batch_6.txt}; do
    sbatch -J "annotate_RefSeq_${list}" \
           -o ${DIR}/logs/annotate_coord_${list}_all.slog \
           -t 00-120 \
           --cpus-per-task 34 \
           --mem-per-cpu=24G \
        --wrap="cat ${list} \
            | /usr/bin/time -v $METAGRAPH annotate -v \
                -i ${DIR}/../../graph_k31.dbg \
                -p 2 \
                --mem-cap-gb 10 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                --threads-each 17 \
                --anno-filename \
                --coordinates \
                --separately \
                -o ${DIR}/columns_new"; \
done


DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
pushd ${DIR}/batches
cp ${DIR}/batches/input/input.txt columns.txt
sed -i 's/raw_data/data/' columns.txt
sed -i 's/fna\/complete_split/output_k31_annotation\/k31_taxid\/nobackup\/build\/columns/' \
    columns.txt
sed -i 's/$/.column.annodbg/' columns.txt
split -n r/6 columns.txt
popd


DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
mkdir ${DIR}/rd;
mkdir ${DIR}/rd/rd_columns;
ln -s ${DIR}/../../graph.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
mkdir ${DIR}/rd/rd_columns;
cd ${DIR}/batches;
for list in x*; do
    bsub -J "RefSeq_rd_0_${list}" \
         -oo ${DIR}/logs/rd_0_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19500] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff --coordinates \
                --max-path-length 100 \
                --row-diff-stage 0 \
                --mem-cap-gb 450 \
                -p 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg"; \
done

for list in x*; do
    bsub -J "RefSeq_rd_1_${list}" \
         -w "RefSeq_rd_0_${list}" \
         -oo ${DIR}/logs/rd_1_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19500] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff --coordinates \
                --max-path-length 100 \
                --row-diff-stage 1 \
                --mem-cap-gb 450 \
                -p 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg"; \
done

DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
cd ${DIR}/batches;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for list in x*; do
    bsub -J "RefSeq_rd_2_${list}" \
         -oo ${DIR}/logs/rd_2_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19500] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff --coordinates \
                --max-path-length 100 \
                --row-diff-stage 2 \
                --mem-cap-gb 450 \
                -p 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg"; \
done


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
bsub -J "refseq_rd_coord" \
     -w "RefSeq_rd_2_*" \
     -oo ${DIR}/logs/rd_coord.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19500] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_coord \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
bsub -J "refseq_rd_brwt_coord" \
     -oo ${DIR}/logs/rd_brwt_coord_build.lsf \
     -W 72:00 \
     -n 48 -R "rusage[mem=20000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt_coord \
            --linkage-file ${DIR}/annotation.linkage \
            --parallel-nodes 12 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 96"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/refseq/release97/output_k31_annotation/k31_taxid/nobackup/build;
bsub -J "refseq_rd_brwt_relax" \
     -oo ${DIR}/logs/build_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19500] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 16 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt_coord.annodbg"

######################## TARA OCEANS ##########################

```bash
DATA=/cluster/work/grlab/projects/metagenome/raw_data/tara/assemblies-fasta
DIR=~/metagenome/data/tara/assemblies;

mkdir $DIR;
mkdir $DIR/logs;

find $DATA -name "*.gz" > $DIR/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_graph_assemblies" \
     -oo $DIR/logs/build_graph.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat $DIR/list.txt \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/graph; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 72 \
            -o $DIR/graph_small \
            $DIR/graph.dbg;"

DATA=/cluster/work/grlab/projects/metagenome/raw_data/tara/assemblies-fasta
DIR=~/metagenome/data/tara/assemblies;
mkdir $DIR/columns;
mkdir $DIR/batches;
cd $DIR/batches;
split -d -n r/40 <(find $DATA -name "*.gz" | shuf);

bsub -J "count_bp_assemblies[1-$(cat $DIR/list.txt | wc -l)]" \
     -o /dev/null \
     -W 4:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
    "file=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/list.txt); \
    id=\\\$(basename \\\${file%.fasta.gz}); \
    echo \\\${id} \\\$(zcat \\\$file | sed '/^>/d' | tr -d '\n' | wc -c) >> $DIR/num_bp.txt;";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "annotate_assemblies_${list}" \
         -w "build_graph_assemblies" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=15000] span[hosts=1]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/graph.dbg \
                --anno-header \
                --separately \
                --coordinates \
                -o ${DIR}/columns \
                -p 4 \
                --threads-each 9"; \
done

DIR=~/metagenome/data/tara/assemblies;
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/graph.dbg ${DIR}/rd/graph.dbg;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/tara/assemblies;
sbatch -J "assemblies_rd_0" \
       -o ${DIR}/logs/rd_0.slog \
       -t 00-72 \
       --cpus-per-task 34 \
       --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 68";

DIR=~/metagenome/data/tara/assemblies;
sbatch -J "assemblies_rd_1" \
       -d afterok:$(get_jobid assemblies_rd_0) \
       -o ${DIR}/logs/rd_1.slog \
       -t 00-72 \
       --cpus-per-task 34 \
       --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 300 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 68";

DIR=~/metagenome/data/tara/assemblies;
sbatch -J "assemblies_rd_2" \
       -d afterok:$(get_jobid assemblies_rd_1) \
       -o ${DIR}/logs/rd_2.slog \
       -t 00-24 \
       --cpus-per-task 34 \
       --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 300 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 68";

sbatch -J "assemblies_rd_disk" \
     -d afterok:$(get_jobid assemblies_rd_2) \
     -o $DIR/logs/rd_disk.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --partition=compute,gpu \
     --wrap="find $DIR/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 220 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

DIR=~/metagenome/data/tara/assemblies;
sbatch -J "assemblies_rd_sparse" \
       -d afterok:$(get_jobid assemblies_rd_2) \
       -o ${DIR}/logs/rd_sparse.slog \
       -t 00-24 \
       --cpus-per-task 36 \
       --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72";

DIR=~/metagenome/data/tara/assemblies;
sbatch -J "assemblies_rd_flat" \
       -d afterok:$(get_jobid assemblies_rd_2) \
       -o ${DIR}/logs/rd_flat.slog \
       -t 00-24 \
       --cpus-per-task 34 \
       --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 68";

DIR=~/metagenome/data/tara/assemblies;
bsub -J "assemblies_rd_brwt" \
     -w "assemblies_rd_2" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10";

DIR=~/metagenome/data/tara/assemblies;
bsub -J "assemblies_rd_brwt_relax" \
     -w "assemblies_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=10000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt.annodbg";
```


######################## TARA OCEANS (genome files) ##########################

```bash
DATA=/cluster/work/grlab/projects/metagenome/raw_data/tara/genomes-fasta;
DIR=~/metagenome/data/tara/genome_coord;

mkdir $DIR;
mkdir $DIR/logs;

find $DATA -name "*.gz" > $DIR/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_graph_genomes" \
     -oo $DIR/logs/build_graph.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat $DIR/list.txt \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/graph; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 72 \
            -o $DIR/graph_small \
            $DIR/graph.dbg;"

DATA=/cluster/work/grlab/projects/metagenome/raw_data/tara/genomes-fasta
DIR=~/metagenome/data/tara/genome_coord;
mkdir $DIR/columns;
mkdir $DIR/batches;
cd $DIR/batches;
split -d -n r/40 <(find $DATA -name "*.gz" | shuf);

for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "count_bp_assemblies_${list}" \
         -o /dev/null \
         -W 4:00 \
         -n 1 -R "rusage[mem=10000] span[hosts=1]" \
    "for file in \\\$(cat $DIR/batches/${list}); do \
        id=\\\$(basename \\\${file%.fa.gz}); \
        echo \\\${id} \\\$(zcat \\\$file | sed '/^>/d' | tr -d '\n' | wc -c) >> $DIR/num_bp.txt; \
    done";
done

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "annotate_genomes_${list}" \
         -w "build_graph_genomes" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=15000] span[hosts=1]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/graph.dbg \
                --anno-filename \
                --separately \
                --coordinates \
                -o ${DIR}/columns \
                -p 4 \
                --threads-each 9"; \
done

DIR=~/metagenome/data/tara/genome_coord;
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/graph.dbg ${DIR}/rd/graph.dbg;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/tara/genome_coord;
bsub -J "genomes_rd_0" \
     -w "annotate_genomes_*" \
     -o ${DIR}/logs/rd_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --coordinates \
            --mem-cap-gb 600 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

DIR=~/metagenome/data/tara/genome_coord;
bsub -J "genomes_rd_1" \
     -w "genomes_rd_0" \
     -o ${DIR}/logs/rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --coordinates \
            --mem-cap-gb 300 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

DIR=~/metagenome/data/tara/genome_coord;
bsub -J "genomes_rd_2" \
     -w "genomes_rd_1" \
     -oo ${DIR}/logs/rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --coordinates \
            --mem-cap-gb 300 \
            --disk-swap \\\"\\\" \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

DIR=~/metagenome/data/tara/genome_coord;
bsub -J "genomes_rd_brwt_coord" \
     -w "genomes_rd_2" \
     -oo ${DIR}/logs/rd_brwt_coord.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt_coord \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10";

DIR=~/metagenome/data/tara/genome_coord;
bsub -J "genomes_rd_brwt_coord_relax" \
     -w "genomes_rd_brwt_coord" \
     -oo ${DIR}/logs/rd_brwt_coord_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=10000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt_coord.annodbg";
```


######################## MetaGut ##########################

rm -rf ~/metagenome/data/cloudcompute/metagut_graphs/build/rd;
mkdir ~/metagenome/data/cloudcompute/metagut_graphs/build/rd;

ln -s ~/metagenome/data/cloudcompute/metagut_graphs/graph_merged_complete_k31.primary.dbg \
      ~/metagenome/data/cloudcompute/metagut_graphs/build/rd/graph.dbg;

cp ~/metagenome/data/cloudcompute/metagut_annotation/files_to_annotate.txt \
   ~/metagenome/data/cloudcompute/metagut_graphs/build/columns.txt
sed -i 's/^.*\///' ~/metagenome/data/cloudcompute/metagut_graphs/build/columns.txt
sed -i 's/^/\/cluster\/home\/mikhaika\/metagenome\/data\/cloudcompute\/metagut_annotation\/columns\//' \
            ~/metagenome/data/cloudcompute/metagut_graphs/build/columns.txt
sed -i 's/$/.column.annodbg/' ~/metagenome/data/cloudcompute/metagut_graphs/build/columns.txt

DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
mkdir ${DIR}/batches;
pushd ${DIR}/batches
cp ../columns.txt ./
split -n r/6 columns.txt
popd

DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
mkdir ${DIR}/rd/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
cd ${DIR}/batches;
for list in xa*; do
    bsub -J "MetaGut_rd_stage_0_${list}" \
         -oo ${DIR}/logs/rd_stage_0_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 0 \
                --mem-cap-gb 650 \
                --parallel 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg \
                2>&1 | tee ${DIR}/logs/rd_stage_0_${list}.log"; \
done

bsub -J "MetaGut_rd_stage_1" \
     -oo ${DIR}/logs/rd_stage_1.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat /dev/null \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 650 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_1.log"; \

for list in xa*; do
    bsub -J "MetaGut_rd_stage_1_${list}" \
         -oo ${DIR}/logs/rd_stage_1_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 1 \
                --mem-cap-gb 500 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                --parallel 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg \
                2>&1 | tee ${DIR}/logs/rd_stage_1_${list}.log"; \
done

bsub -J "MetaGut_rd_stage_2" \
     -oo ${DIR}/logs/rd_stage_2.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat /dev/null \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 100 \
            --mem-cap-gb 650 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_2.log"; \

for list in xa*; do
    bsub -J "MetaGut_rd_stage_2_${list}" \
         -w "MetaGut_rd_stage_2" \
         -oo ${DIR}/logs/rd_stage_2_${list}.lsf \
         -W 72:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${list} \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 2 \
                --max-path-length 100 \
                --mem-cap-gb 650 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                --parallel 72 \
                -o ${DIR}/rd/rd_columns/out \
                -i ${DIR}/rd/graph.dbg \
                2>&1 | tee ${DIR}/logs/rd_stage_2_${list}.log"; \
done


DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
bsub -J "metagut_cluster" \
     -oo ${DIR}/logs/cluster_original.lsf \
     -W 240:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat ${DIR}/columns.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type brwt \
            --linkage --greedy --subsample 50000 \
            -o ${DIR}/cluster_original \
            -p 72 \
            2>&1 | tee ${DIR}/logs/cluster_original.log"

DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
for x in $(cat ${DIR}/columns.txt); do
    x=$(basename $x);
    echo $DIR/rd/rd_columns/${x%.column.annodbg}.row_diff.annodbg;
done > ${DIR}/rd_columns_cluster_original.txt
bsub -J "metagut_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 48:00 \
     -n 48 -R "rusage[mem=40400] span[hosts=1]" \
    "cat ${DIR}/rd_columns_cluster_original.txt \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --linkage-file ${DIR}/cluster_original \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 96 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/rd_brwt.log"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
bsub -J "metagut_rd_brwt_relax" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 48:00 \
     -n 18 -R "rusage[mem=78000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/rd_brwt_relax.log"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
sbatch -J "metagut_rd_disk" \
     -o $DIR/logs/rd_disk.slog \
     -t 00-120 \
     --cpus-per-task 56 \
     --mem-per-cpu=15G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 800 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 56"

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/cloudcompute/metagut_graphs/nobackup/build;
sbatch -J "metagut_rd_disk_big" \
     -o $DIR/logs/rd_disk_big.slog \
     -t 00-120 \
     --cpus-per-task 70 \
     --mem-per-cpu=20G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 1200 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation_big \
            -p 70"



######################## UHGG small ##########################

ln -s ~/metagenome/data/uhgg/uhgg_catalogue/graphs/graph_complete_k31.dbg \
      ~/metagenome/data/uhgg/uhgg_catalogue/graphs/build/rd/graph.dbg;

DIR=~/metagenome/data/uhgg/uhgg_catalogue/graphs/build;
mkdir ${DIR}/rd/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

bsub -J "uhgg_rd_stage_0" \
     -oo ${DIR}/logs/rd_stage_0.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=5000] span[hosts=1]" \
    "find ~/metagenome/data/uhgg/uhgg_catalogue/annotation/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 650 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_0.log"; \

bsub -J "uhgg_rd_stage_1" \
     -w "uhgg_rd_stage_0" \
     -oo ${DIR}/logs/rd_stage_1.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=5000] span[hosts=1]" \
    "find ~/metagenome/data/uhgg/uhgg_catalogue/annotation/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 650 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_1.log"; \

bsub -J "uhgg_rd_stage_2" \
     -w "uhgg_rd_stage_1" \
     -oo ${DIR}/logs/rd_stage_2.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=5000] span[hosts=1]" \
    "find ~/metagenome/data/uhgg/uhgg_catalogue/annotation/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 200 \
            --mem-cap-gb 650 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_2.log"; \

bsub -J "uhgg_rd_brwt" \
     -w "uhgg_rd_stage_2" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=5000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 10000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/rd_brwt.log"

bsub -J "uhgg_rd_brwt_relax" \
     -w "uhgg_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=5000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        -p 72 \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/rd_brwt_relax.log"

######################## UHGG ##########################

rm -rf ~/metagenome/data/uhgg/all_genomes/build/rd;
mkdir ~/metagenome/data/uhgg/all_genomes/build/rd;

ln -s ~/metagenome/data/uhgg/all_genomes/graph_complete_k31.dbg \
      ~/metagenome/data/uhgg/all_genomes/build/rd/graph.dbg;

DIR=~/metagenome/data/uhgg/all_genomes/build;
mkdir ${DIR}/rd/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

bsub -J "UHGG_rd_stage_0" \
     -oo ${DIR}/logs/rd_stage_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 650 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_0.log"; \

bsub -J "UHGG_rd_stage_1" \
     -w "UHGG_rd_stage_0" \
     -oo ${DIR}/logs/rd_stage_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 650 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_1.log"; \

bsub -J "UHGG_rd_stage_2" \
     -w "UHGG_rd_stage_1" \
     -oo ${DIR}/logs/rd_stage_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 200 \
            --mem-cap-gb 650 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 72 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg \
            2>&1 | tee ${DIR}/logs/rd_stage_2.log"; \

sbatch -J "UHGG_rd_disk" \
     -d afterok:$(get_jobid UHGG_rd_stage_2) \
     -o $DIR/logs/rd_disk.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
     --partition=compute,gpu \
     --wrap="find $DIR/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 220 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

bsub -J "UHGG_rd_brwt" \
     -w "UHGG_rd_stage_2" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/rd_brwt.log"

bsub -J "UHGG_rd_brwt_relax" \
     -w "UHGG_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        -p 72 \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/rd_brwt_relax.log"


######################## MetaSUB 19 ##########################

mkdir ~/metagenome/data/metasub/graphs/k19/build/rd;

ln -s ~/metagenome/data/metasub/graphs/k19/graph_merged_k19.primary.dbg \
      ~/metagenome/data/metasub/graphs/k19/build/rd/graph.dbg;

DIR=~/metagenome/data/metasub/graphs/k19/build;
mkdir ${DIR}/rd/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

sbatch -J "metasub_rd_0" \
     -o $DIR/logs/rd_0.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=14G \
     --partition=compute \
    --wrap="find $DIR/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 420 \
            -p 34 \
            -o $DIR/rd/rd_columns/out \
            -i $DIR/rd/graph.dbg"; \

sbatch -J "metasub_rd_1" \
     -d afterok:$(get_jobid metasub_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=14G \
     --partition=compute \
    --wrap="find $DIR/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 420 \
            -p 34 \
            -o $DIR/rd/rd_columns/out \
            -i $DIR/rd/graph.dbg"; \

sbatch -J "metasub_rd_2" \
     -d afterok:$(get_jobid metasub_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=14G \
     --partition=compute \
    --wrap="find $DIR/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 420 \
            -p 34 \
            -o $DIR/rd/rd_columns/out \
            -i $DIR/rd/graph.dbg"; \

sbatch -J "metasub_rd_disk" \
     -d afterok:$(get_jobid metasub_rd_2) \
     -o $DIR/logs/rd_disk.slog \
     -t 0-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --wrap="find $DIR/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 420 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

bsub -J "metasub_rd_brwt" \
     -w "metasub_rd_stage_2" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=15000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 10000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/rd_brwt.log"

bsub -J "metasub_rd_brwt_relax" \
     -w "metasub_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=15000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        -p 72 \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/rd_brwt_relax.log"


######################## MetaSUB 41 ##########################

ln -s ~/metagenome/data/metasub/graphs/output_k41_cleaned_graph/graph_merged_k41.primary.dbg \
      ~/metagenome/data/metasub/graphs/k41/build/rd/graph.dbg;

DIR=~/metagenome/data/metasub/graphs/k41/build;
mkdir ${DIR}/rd/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

sbatch -J "metasub_rd_stage_0" \
     -o ${DIR}/logs/rd_stage_0.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=19G \
     --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 300 \
            --parallel 36 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg"; \

sbatch -J "metasub_rd_stage_1" \
     -d afterok:$(get_jobid metasub_rd_stage_0) \
     -o ${DIR}/logs/rd_stage_1.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=19G \
     --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 300 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 36 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg"; \

sbatch -J "metasub_rd_stage_2" \
     -d afterok:$(get_jobid metasub_rd_stage_1) \
     -o ${DIR}/logs/rd_stage_2.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=19G \
     --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 100 \
            --mem-cap-gb 300 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --parallel 36 \
            -o ${DIR}/rd/rd_columns/out \
            -i ${DIR}/rd/graph.dbg"; \

sbatch -J "metasub_rd_brwt" \
     -d afterok:$(get_jobid metasub_rd_stage_2) \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=15G \
     --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 10000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10"

sbatch -J "metasub_rd_brwt_relax" \
     -d afterok:$(get_jobid metasub_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=15G \
     --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        -p 72 \
        ${DIR}/annotation.row_diff_brwt.annodbg"

sbatch -J "metasub_rd_disk" \
     -d afterok:$(get_jobid metasub_rd_stage_2) \
     -o ${DIR}/logs/rd_disk.slog \
     -t 00-24 \
     --cpus-per-task 36 \
     --mem-per-cpu=15G \
     --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 300 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 36"


######################## FUNGI ##########################
```bash
DIR=~/metagenome/data/cloudcompute/fungi_graphs;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

mkdir -p ${DIR}/logs;

sbatch -J "fungi" \
     -o $DIR/logs/build_graph.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/samples.txt \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode canonical \
            --mem-cap-gb 80 \
            --disk-swap ~/metagenome/scratch/nobackup/ \
            -p 34 \
            -o $DIR/graph; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta \
            --primary-kmers \
            -p 34 \
            -o $DIR/primary_contigs \
            $DIR/graph.dbg; \
    /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode primary \
            --mem-cap-gb 80 \
            --disk-swap ~/metagenome/scratch/nobackup/ \
            -p 34 \
            -o $DIR/graph_primary \
            $DIR/primary_contigs.fasta.gz; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 34 \
            -o $DIR/graph_primary_small \
            $DIR/graph_primary.dbg;"

mkdir ${DIR}/columns;

sbatch -J "fungi_annotate" \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/samples.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 4 \
            --threads-each 8"

mkdir -p ${DIR}/rd/rd_columns;

ln -s $DIR/graph_primary.dbg ${DIR}/rd/graph.dbg;

sbatch -J "fungi_rd_0" \
     -o $DIR/logs/rd_0.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --exclude compute-biomed-03,compute-biomed-02,compute-biomed-10 \
    --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            -p 34 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out"

sbatch -J "fungi_rd_1" \
     -d afterok:$(get_jobid fungi_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --exclude compute-biomed-03,compute-biomed-02,compute-biomed-10 \
    --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out"

sbatch -J "fungi_rd_2" \
     -d afterok:$(get_jobid fungi_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
     --exclude compute-biomed-03,compute-biomed-02,compute-biomed-10 \
    --wrap="find ${DIR}/columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out"

sbatch -J "fungi_rd_flat" \
     -d afterok:$(get_jobid fungi_rd_2) \
     -o $DIR/logs/rd_flat.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
     --exclude compute-biomed-03,compute-biomed-02 \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "fungi_rd_sparse" \
     -d afterok:$(get_jobid fungi_rd_2) \
     -o $DIR/logs/rd_sparse.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
     --exclude compute-biomed-03,compute-biomed-02 \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "fungi_rd_disk" \
     -d afterok:$(get_jobid fungi_rd_2) \
     -o $DIR/logs/rd_disk.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=15G \
     --exclude compute-biomed-03,compute-biomed-02 \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 300 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "fungi_rd_brwt" \
     -d afterok:$(get_jobid fungi_rd_2) \
     -o $DIR/logs/rd_brwt.slog \
     -t 00-72 \
     --cpus-per-task 56 \
     --mem-per-cpu=15G \
     --exclude compute-biomed-03,compute-biomed-02 \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 200000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation_new \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 56 --parallel-nodes 10"

sbatch -J "fungi_rd_brwt_relax" \
     -d afterok:$(get_jobid fungi_rd_brwt) \
     -o $DIR/logs/rd_brwt_relax.slog \
     -t 00-72 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
     --exclude compute-biomed-03,compute-biomed-02 \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation_new.relaxed \
        ${DIR}/annotation_new.row_diff_brwt.annodbg"
```

######################## FUNGI SUBSET ##########################

DIR=~/metagenome/data/cloudcompute/fungi5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir -p ${DIR}/logs;

sbatch -J "fungi5K" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode canonical \
            -p 34 \
            -o $DIR/graph; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta \
            --primary-kmers \
            -p 34 \
            -o $DIR/primary_contigs \
            $DIR/graph.dbg; \
    /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode primary \
            -p 34 \
            -o $DIR/graph_primary \
            $DIR/primary_contigs.fasta.gz; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 34 \
            -o $DIR/graph_primary_small \
            $DIR/graph_primary.dbg;"


mkdir ${DIR}/columns;

sbatch -J "fungi5K_annotate" \
     -d afterok:$(get_jobid fungi5K) \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 8 \
            --threads-each 8"


rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "fungi5K_rd_0" \
     -d afterok:$(get_jobid fungi5K_annotate) \
     -o $DIR/logs/rd_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 150 \
            -p 34"


sbatch -J "fungi5K_rd_1" \
     -d afterok:$(get_jobid fungi5K_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 150 \
            -p 34"


sbatch -J "fungi5K_rd_2" \
     -d afterok:$(get_jobid fungi5K_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 150 \
            -p 34"

DIR=~/metagenome/data/cloudcompute/fungi5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

sbatch -J "fungi5K_rd_flat" \
     -d afterok:$(get_jobid fungi5K_rd_2) \
     -o ${DIR}/logs/rd_flat.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "fungi5K_rd_sparse" \
     -d afterok:$(get_jobid fungi5K_rd_2) \
     -o ${DIR}/logs/rd_sparse.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "fungi5K_rd_brwt" \
     -d afterok:$(get_jobid fungi5K_rd_2) \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 10"

sbatch -J "fungi5K_rd_brwt_relax" \
     -d afterok:$(get_jobid fungi5K_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"


######################## PLANTS ##########################

DIR=~/metagenome/data/cloudcompute/viridiplantae_graphs;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
bsub -J "plants_primary" \
     -oo ${DIR}/build_primary_new.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonE7_8867v3]" \
    "gtime -v $METAGRAPH build -v \
            -k 31 \
            --mode primary \
            -o $DIR/graph_primary_new \
            $DIR/primary_contigs_k31.fasta.gz \
            --disk-swap ~/metagenome/scratch/nobackup/ \
            --mem-cap-gb 100 \
            --disk-cap-gb 20000 \
            -p 36 \
            2>&1 | tee ${DIR}/build_primary_new.log";


ln -s ~/metagenome/data/cloudcompute/viridiplantae_graphs/graph_primary.dbg \
      ~/metagenome/data/cloudcompute/viridiplantae_annotation/rd/graph_primary.dbg
mkdir ~/metagenome/data/cloudcompute/viridiplantae_annotation/rd/rd_columns

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation;
cd ${DIR}/rd/batches;
for S in {xac,}; do
    bsub -J "plants_rd_${S}" \
         -oo ${DIR}/rd/transform_to_rd_${S}.lsf \
         -W 120:00 \
         -n 36 -R "rusage[mem=20000] span[hosts=1] select[model==XeonE7_8867v3]" \
        "cat ${DIR}/rd/batches/${S} \
            | gtime -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length 100 \
                --mem-cap-gb 650 \
                -i ${DIR}/rd/graph_primary.dbg \
                -o ${DIR}/rd/rd_columns/out \
                -p 72 \
                2>&1 | tee ${DIR}/rd/transform_to_rd_${S}.log";
done


cd ~/metagenome/data/cloudcompute/viridiplantae_annotation/rd/second_run;
sed -i 's/.fasta.gz.column.annodbg//' all_columns.txt
sed -i 's/.fasta.gz.row_diff.annodbg//' rd_transformed.txt
comm -23 all_columns.txt rd_transformed.txt > remaining.txt

sed -i 's/$/.fasta.gz.column.annodbg/' remaining.txt
sed -i 's/^/\/cluster\/home\/mikhaika\/metagenome\/data\/cloudcompute\/viridiplantae_annotation\/columns\//' remaining.txt
cat remaining.txt | shuf > remaining_random.txt
split -n r/6 remaining_random.txt


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation;
cd ${DIR}/rd/second_run;
for S in xa*; do
    bsub -J "plants_rd_${S}" \
         -oo ${DIR}/rd/transform_to_rd_run2_${S}.lsf \
         -W 120:00 \
         -n 36 -R "rusage[mem=19500] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${DIR}/rd/second_run/${S} \
            | gtime -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length 100 \
                --mem-cap-gb 650 \
                -i ${DIR}/rd/graph_primary.dbg \
                -o ${DIR}/rd/rd_columns/out \
                -p 72 \
                2>&1 | tee ${DIR}/rd/transform_to_rd_run2_${S}.log";
done


sed -i 's/$/.fasta.gz.column.annodbg/' all_columns.txt
sed -i 's/^/\/cluster\/home\/mikhaika\/metagenome\/data\/cloudcompute\/viridiplantae_annotation\/columns\//' all_columns.txt
split -n r/12 <(cat all_columns.txt | shuf)

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation/nobackup/build;
mkdir ${DIR}/rd/rd_columns_opt;
cd ${DIR}/rd/rd_columns_opt;
ln -s ../rd_columns/*row_reduction* .;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation/nobackup/build;
cd ${DIR}/batches/left;
for S in x*; do
    bsub -J "plants_rd_${S}_opt" \
         -oo ${DIR}/logs/transform_to_rd_${S}__opt.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=19500] span[hosts=1] select[hname!='le-amd-fp-004']" \
        "cat ${DIR}/batches/left/${S} \
            | gtime -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length 100 \
                --mem-cap-gb 650 \
                --optimize \
                -i ${DIR}/rd/graph_primary.dbg \
                -o ${DIR}/rd/rd_columns_opt/out \
                -p 72 \
                2>&1 | tee ${DIR}/logs/transform_to_rd_${S}__opt.log";
done

cp columns.viridiplantae.filtered.txt rd_columns.txt
sed -i 's/$/.fasta.gz.row_diff.annodbg/' rd_columns.txt
sed -i 's/^/\/cluster\/work\/grlab\/projects\/metagenome\/data\/cloudcompute\/viridiplantae_annotation\/nobackup\/build\/rd\/rd_columns_opt\//' \
    rd_columns.txt


bsub -J "plants_rd_brwt_reserve" \
     -oo /dev/null \
     -W 120:00 \
     -n 16 -R "rusage[mem=4000] span[hosts=1] select[hname=='le-fat-001']" \
    "sleep 1000h";
bsub -J "plants_rd_brwt_reserve" \
     -oo /dev/null \
     -W 120:00 \
     -n 48 -R "rusage[mem=4000] span[hosts=1] select[hname=='le-fat-001']" \
    "sleep 1000h";


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation/nobackup/build;
bsub -J "plants_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt.lsf \
     -W 48:00 \
     -n 48 -R "rusage[mem=40000] span[hosts=1]" \
    "cat ${DIR}/rd_columns.txt \
        | gtime -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            -i ${DIR}/rd/graph_primary.dbg \
            --linkage-file ${DIR}/columns.viridiplantae.filtered.taxids.linkage.txt \
            -o ${DIR}/annotation \
            -p 96 --parallel-nodes 5 \
            2>&1 | tee ${DIR}/logs/build_rd_brwt.log"


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/cloudcompute/viridiplantae_annotation/nobackup/build;
bsub -J "plants_rd_brwt_relax" \
     -w "plants_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt_relax.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=15000] span[hosts=1]" \
    "gtime -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 32 \
        -o ${DIR}/../../annotation.relaxed \
        ${DIR}/../../annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/build_rd_brwt_relax.log"

```


```bash
graph=/cluster/work/grlab/projects/metagenome/data/uhgg/uhgg_catalogue/graphs/graph_complete_k31.dbg
anno=/cluster/work/grlab/projects/metagenome/data/uhgg/uhgg_catalogue/annotation/output_k31.relaxed.relabeled.row_diff_brwt.annodbg
out_dir=${graph}.benchmarks
mkdir $out_dir
for QUERY in ~/metagenome/data/BIGSI/subsets/query/samples/haib18CEM5453_HMCMJCCXY_SL336225.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/nucleotide_fasta_protein_homolog_model.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/DRR067889.fasta; do
    /cluster/home/mikhaika/projects/projects2014-metagenome/metagraph/scripts/run_query.sh "$graph" "$anno" "$QUERY" "$out_dir";
done
```


```bash
######################## HUMAN ##########################

DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;

rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s ~/metagenome/data/cloudcompute/homo_sapiens_graphs/graph_merged_complete_k31.primary.dbg \
      ${DIR}/rd/graph.dbg;

cp ~/metagenome/data/cloudcompute/homo_sapiens_annotation/files_to_annotate.txt \
   ${DIR}/columns.txt
sed -i 's/^.*\///' ${DIR}/columns.txt
sed -i 's/^/\/cluster\/home\/mikhaika\/metagenome\/data\/cloudcompute\/homo_sapiens_annotation\/columns\//' \
            ${DIR}/columns.txt
sed -i 's/$/.column.annodbg/' ${DIR}/columns.txt


mkdir ${DIR}/batches;
pushd ${DIR}/batches
cp ../columns.txt ./
split -l -n r/60 columns.txt
popd

mkdir ${DIR}/rd/rd_columns

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph

for B in {64..79..4}; do
    dep_cond="";
    for N in `seq $B $((B+3))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "human_rd_stage_0_$CHUNK" \
                 -oo ${DIR}/logs/human_rd_stage_0_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 48:00 \
                 -n 36 -R "rusage[mem=13000] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 0 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_count \
                    --mem-cap-gb 450 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/human_rd_stage_0_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done

bsub -J "human_rd_stage_1" \
     -oo ${DIR}/logs/human_rd_stage_1.lsf \
     -w "1322561 && 1322569 && 1322577 && 1322585 && 1322593 && 1322601 && 1322609 && 1322617 && 1322641 && 1322645 && 1322649 && 1322653 && 1322562 && 1322563 && 1322564 && 1322565 && 1322566 && 1322567 && 1322568 && 1322570 && 1322571 && 1322572 && 1322573 && 1322574 && 1322575 && 1322576 && 1322578 && 1322579 && 1322580 && 1322581 && 1322582 && 1322583 && 1322584 && 1322586 && 1322587 && 1322588 && 1322589 && 1322590 && 1322591 && 1322592 && 1322594 && 1322595 && 1322596 && 1322597 && 1322598 && 1322599 && 1322600 && 1322602 && 1322603 && 1322604 && 1322605 && 1322606 && 1322607 && 1322608 && 1322610 && 1322611 && 1322612 && 1322613 && 1322614 && 1322615 && 1322616 && 1322618 && 1322619 && 1322620 && 1322621 && 1322622 && 1322623 && 1322624 && 1322642 && 1322643 && 1322644 && 1322646 && 1322647 && 1322648 && 1322650 && 1322651 && 1322652 && 1322654 && 1322655 && 1322656" \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
    "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
        --anno-type row_diff \
        --row-diff-stage 1 \
        -i ${DIR}/rd/graph.dbg \
        -o ${DIR}/rd/rd_columns/vectors \
        --mem-cap-gb 600 \
        -p 72 \
        2>&1 | tee ${DIR}/logs/human_rd_stage_1.log"


DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph
rm ${DIR}/rd/rd_columns/*.row_count

# 4, 48:00, mem=13400, 335
# 8, 48:00, mem=19000, 475
# 8, 48:00, mem=20000, 500
# 6, 48:00, mem=20000, 450

for B in {64..79..4}; do
    dep_cond="";
    for N in `seq $B $((B+3))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "human_rd_stage_1_$CHUNK" \
                 -oo ${DIR}/logs/human_rd_stage_1_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 48:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 1 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_reduction \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --mem-cap-gb 500 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/human_rd_stage_1_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done


DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph
bsub -J "human_rd_stage_2" \
     -oo ${DIR}/logs/human_rd_stage_2.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19400] span[hosts=1]" \
        "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 100 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 650 \
            -p 72 \
            2>&1 | tee ${DIR}/logs/human_rd_stage_2.log"


for N in {0..79}; do
    N=$(printf "%02d" $N);
    CHUNK=x$N;
    bsub -J "human_rd_stage_2_$CHUNK" \
         -oo ${DIR}/logs/human_rd_stage_2_$CHUNK.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=13000] span[hosts=1]" \
        "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --max-path-length 100 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 450 \
            -p 72 \
            2>&1 | tee ${DIR}/logs/human_rd_stage_2_$CHUNK.log";
done

DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
L=~/metagenome/data/taxonomy/linkage/homo_sapiens/redo/split_rd;
for G in {0..10}; do
    cp $L/group_$G $DIR/;
    sed -i 's/^/\/cluster\/work\/grlab\/projects\/metagenome\/data\/cloudcompute\/homo_sapiens_graphs\/build\/nobackup\/rd\/rd_columns\//' $DIR/group_$G;
    sed -i 's/$/.fasta.gz.row_diff.annodbg/' $DIR/group_$G;
done

DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
L=~/metagenome/data/taxonomy/linkage/homo_sapiens/redo/split_rd;
for G in {0,1}; do
    cp $L/group_$G.linkage ${DIR}/;
    bsub -J "homo_sapiens_rd_brwt_$G" \
         -oo ${DIR}/logs/rd_brwt_$G.lsf \
         -W 24:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "cat $DIR/group_$G \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_brwt \
                --linkage-file $DIR/group_$G.linkage \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/annotation_$G \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -p 72 --parallel-nodes 10 \
                2>&1 | tee ${DIR}/logs/rd_brwt_$G.log"
done

DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
for G in {5..10}; do
    bsub -J "homo_sapiens_rd_brwt_$G" \
         -oo ${DIR}/logs/rd_brwt_$G.lsf \
         -W 48:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "cat $DIR/group_$G \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_brwt \
                --greedy --subsample 10000000 \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/annotation_$G \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -p 72 --parallel-nodes 10 \
                2>&1 | tee ${DIR}/logs/rd_brwt_$G.log"
done

#################TEST
DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
L=~/metagenome/data/taxonomy/linkage/homo_sapiens/redo/split_rd;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;
G=9
cp $L/group_$G $DIR/group_${G}_test;
sed -i 's/^/\/cluster\/work\/grlab\/projects\/metagenome\/data\/cloudcompute\/homo_sapiens_annotation\/nobackup\/columns\//' $DIR/group_${G}_test;
sed -i 's/$/.fasta.gz.column.annodbg/' $DIR/group_${G}_test;

for G in {9,}; do
    bsub -J "homo_sapiens_rd_brwt_$G" \
         -oo ${DIR}/logs/brwt_${G}.lsf \
         -W 48:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "cat $DIR/group_${G}_test \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type brwt \
                --greedy --subsample 10000000 \
                -o ${DIR}/annotation_${G} \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -p 72 --parallel-nodes 10 \
                2>&1 | tee ${DIR}/logs/brwt_${G}.log"
done
#################TEST

DIR=~/metagenome/data/cloudcompute/homo_sapiens_graphs/build/nobackup;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
for G in {0..10}; do
    bsub -J "homo_sapiens_rd_brwt_${G}_relax" \
         -w "homo_sapiens_rd_brwt_${G}" \
         -oo ${DIR}/logs/rd_brwt_${G}_relax.lsf \
         -W 24:00 \
         -n 12 -R "rusage[mem=50000] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation_${G}.relaxed \
            ${DIR}/annotation_${G}.row_diff_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/rd_brwt_${G}_relax.log"
done
```


```bash
DIR=~/metagenome/data/cloudcompute/metazoa;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for G in {0..5}; do
    bsub -J "metazoa_canonical_${G}" \
         -oo ${DIR}/build_canonical_group_${G}.lsf \
         -W 120:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "cat ~/metagenome/data/cloudcompute/metazoa/metazoa2_groups_5tb/group_${G}_files.txt | \
            /usr/bin/time -v $METAGRAPH build -v \
                -k 31 \
                --mode canonical \
                --inplace \
                -o $DIR/graph_canonical_${G} \
                --mem-cap-gb 80 \
                --disk-swap ~/metagenome/scratch/nobackup/ \
                -p 72";
done

DIR=~/metagenome/data/cloudcompute/metazoa;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for G in {0..5}; do
    bsub -J "metazoa_${G}_primary_contigs" \
         -w "metazoa_canonical_${G}" \
         -oo ${DIR}/extract_contigs_group_${G}.lsf \
         -W 120:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH transform -v \
                --to-fasta \
                --primary-kmers \
                -o $DIR/primary_contigs_${G} \
                $DIR/graph_canonical_${G}.dbg \
                -p 72";
done

DIR=~/metagenome/data/cloudcompute/metazoa;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "metazoa_canonical" \
     -oo ${DIR}/build_canonical.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=80000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --inplace \
            -o $DIR/graph_basic \
            --disk-swap ~/metagenome/scratch/nobackup/ \
            --mem-cap-gb 180 \
            -p 72 \
            $DIR/primary_contigs*.gz";
````


```bash
######################## Metazoa (MOUSE) ##########################

DIR=~/metagenome/data/cloudcompute/mus_musculus/nobackup/build;

rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s ~/metagenome/data/cloudcompute/mus_musculus/metazoa_mus_musculus_primary.dbg \
      ${DIR}/rd/graph.dbg;

find ~/metagenome/data/cloudcompute/mus_musculus_annotation/nobackup/columns -name "*.annodbg" \
    > ${DIR}/columns.txt

mkdir ${DIR}/batches;
pushd ${DIR}/batches
split -d -n r/15 <(cat ../columns.txt | shuf)
popd

mkdir ${DIR}/rd/rd_columns

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph

for B in {0..12..3}; do
    dep_cond="";
    for N in `seq $B $((B+2))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mouse_rd_stage_0_$CHUNK" \
                 -oo ${DIR}/logs/mouse_rd_stage_0_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 0 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_count \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mouse_rd_stage_0_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done

bsub -J "mouse_rd_stage_1" \
     -oo ${DIR}/logs/mouse_rd_stage_1.lsf \
     -w "1341106 && 1341109 && 1341112 && 1341115 && 1341118 && 1341107 && 1341108 && 1341110 && 1341111 && 1341113 && 1341114 && 1341116 && 1341117 && 1341119 && 1341120" \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
    "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
        --anno-type row_diff \
        --row-diff-stage 1 \
        -i ${DIR}/rd/graph.dbg \
        -o ${DIR}/rd/rd_columns/vectors \
        -p 72 \
        2>&1 | tee ${DIR}/logs/mouse_rd_stage_1.log"


for B in {0..12..3}; do
    dep_cond="    1341121";
    for N in `seq $B $((B+2))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mouse_rd_stage_1_$CHUNK" \
                 -oo ${DIR}/logs/mouse_rd_stage_1_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 1 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_reduction \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mouse_rd_stage_1_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done


bsub -J "mouse_rd_stage_2" \
     -oo ${DIR}/logs/mouse_rd_stage_2.lsf \
     -w "1341137 && 1341138 && 1341139 && 1341140 && 1341141 && 1341142 && 1341143 && 1341144 && 1341145 && 1341146 && 1341147 && 1341148 && 1341149 && 1341150 && 1341151" \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
    "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
        --anno-type row_diff \
        --row-diff-stage 2 \
        -i ${DIR}/rd/graph.dbg \
        -o ${DIR}/rd/rd_columns/vectors \
        -p 72 \
        2>&1 | tee ${DIR}/logs/mouse_rd_stage_2.log"


for B in {0..12..3}; do
    dep_cond="    1341162";
    for N in `seq $B $((B+2))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mouse_rd_stage_2_$CHUNK" \
                 -oo ${DIR}/logs/mouse_rd_stage_2_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 2 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/columns \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mouse_rd_stage_2_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done


find ${DIR}/rd/rd_columns/ -name "*.annodbg" > ${DIR}/rd/rd_columns.txt

bsub -J "mouse_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/build_rd_brwt.log"

bsub -J "mouse_rd_brwt_relax" \
     -w "mouse_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/build_rd_brwt_relax.log"
```



```bash
######################## Metazoa (10k subset) ##########################

DIR=~/metagenome/data/cloudcompute/metazoa/shuf_splits;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
chunk=xan;

bsub -J "metazoa_${chunk}" \
     -oo $DIR/logs/build_graph_${chunk}.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=10000] span[hosts=1]" \
    "cat $DIR/xam \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode canonical \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 40 \
            -o $DIR/graph_${chunk}; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -p 40 \
            -o $DIR/graph_${chunk}_small \
            $DIR/graph_${chunk}.dbg;"


DIR=~/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_10k;

rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s ~/metagenome/data/cloudcompute/metazoa_graphs/metazoa_all_files.shuf.head10000.primary.dbg \
      ${DIR}/rd/graph.dbg;

find ~/metagenome/data/cloudcompute/metazoa_annotation/nobackup/columns.shuf.head10000 -name "*.annodbg" \
    > ${DIR}/columns.txt

mkdir ${DIR}/batches;
pushd ${DIR}/batches
split -d -n r/4 <(cat ../columns.txt | shuf)
popd

mkdir ${DIR}/rd/rd_columns

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph

for B in {0..3..1}; do
    dep_cond="";
    for N in `seq $B $((B+0))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mzoa_head10k_rd_stage_0_$CHUNK" \
                 -oo ${DIR}/logs/mzoa_head10k_rd_stage_0_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 0 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_count \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mzoa_head10k_rd_stage_0_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done

bsub -J "mzoa_head10k_rd_stage_1" \
     -oo ${DIR}/logs/mzoa_head10k_rd_stage_1.lsf \
     -w "1398945 && 1398946 && 1398947 && 1398948" \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
    "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
        --anno-type row_diff \
        --row-diff-stage 1 \
        -i ${DIR}/rd/graph.dbg \
        -o ${DIR}/rd/rd_columns/vectors \
        -p 72 \
        2>&1 | tee ${DIR}/logs/mzoa_head10k_rd_stage_1.log"


for B in {0..3..1}; do
    dep_cond="    1398949";
    for N in `seq $B $((B+0))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mzoa_head10k_rd_stage_1_$CHUNK" \
                 -oo ${DIR}/logs/mzoa_head10k_rd_stage_1_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 1 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/vector_${B}.row_reduction \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mzoa_head10k_rd_stage_1_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done


bsub -J "mzoa_head10k_rd_stage_2" \
     -oo ${DIR}/logs/mzoa_head10k_rd_stage_2.lsf \
     -w " 1398950 && 1398951 && 1398952 && 1398953" \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[hname!='le-amd-fp-004']" \
    "cat /dev/null | /usr/bin/time -v $METAGRAPH transform_anno -v \
        --anno-type row_diff \
        --row-diff-stage 2 \
        -i ${DIR}/rd/graph.dbg \
        -o ${DIR}/rd/rd_columns/vectors \
        -p 72 \
        2>&1 | tee ${DIR}/logs/mzoa_head10k_rd_stage_2.log"


for B in {0..3..1}; do
    dep_cond="    1398954";
    for N in `seq $B $((B+0))`; do
        N=$(printf "%02d" $N);
        CHUNK=x$N;
        JOBID=$(bsub -J "mzoa_head10k_rd_stage_2_$CHUNK" \
                 -oo ${DIR}/logs/mzoa_head10k_rd_stage_2_$CHUNK.lsf \
                 -w "${dep_cond:4}" \
                 -W 24:00 \
                 -n 36 -R "rusage[mem=19400] span[hosts=1] select[hname!='le-amd-fp-004']" \
                "cat ${DIR}/batches/$CHUNK | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --row-diff-stage 2 \
                    -i ${DIR}/rd/graph.dbg \
                    -o ${DIR}/rd/rd_columns/columns \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --mem-cap-gb 650 \
                    -p 72 \
                    2>&1 | tee ${DIR}/logs/mzoa_head10k_rd_stage_2_$CHUNK.log" \
            | awk '/is submitted/{print substr($2, 2, length($2)-2);}');
        dep_cond+=" && $JOBID";
    done
done


find ${DIR}/rd/rd_columns/ -name "*.annodbg" > ${DIR}/rd/rd_columns.txt

bsub -J "mzoa_head10k_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/build_rd_brwt.log"

bsub -J "mzoa_head10k_rd_brwt_relax" \
     -w "mzoa_head10k_rd_brwt" \
     -oo ${DIR}/logs/build_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 36 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg \
        2>&1 | tee ${DIR}/logs/build_rd_brwt_relax.log"
```




```bash
######################## Metazoa (1k studies) ##########################

DIR=~/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_1k_studies;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

mkdir -p ${DIR}/logs;

sbatch -J "metazoa1k" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="cat ${DIR}/samples.txt \
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
            $DIR/graph.dbg; \
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


mkdir ${DIR}/columns;

sbatch -J "metazoa1k_annotate" \
     -d afterok:$(get_jobid metazoa1k) \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/samples.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 5 \
            --threads-each 8"

sbatch -J "metazoa1k_annotate_mid" \
     -o $DIR/logs/annotate_graph_mid.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/samples.txt | tail -n 15000 \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 5 \
            --threads-each 8"

sbatch -J "metazoa1k_annotate_tail" \
     -o $DIR/logs/annotate_graph_tail.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/samples.txt | tail -n 30000 | head -n 15000 \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 5 \
            --threads-each 8"


rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "metazoa1k_rd_0" \
     -d afterok:$(get_jobid metazoa1k_annotate) \
     -o $DIR/logs/rd_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "metazoa1k_rd_1" \
     -d afterok:$(get_jobid metazoa1k_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 650 \
            -p 34 && \
        rm ${DIR}/rd/rd_columns/*.row_count"

sbatch -J "metazoa1k_rd_2" \
     -d afterok:$(get_jobid metazoa1k_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 650 \
            -p 34 && \
        rm ${DIR}/rd/rd_columns/*.row_reduction && \
        rm -r ${DIR}/columns"


sbatch -J "metazoa1k_rd_flat" \
     -d afterok:$(get_jobid metazoa1k_rd_2) \
     -o ${DIR}/logs/rd_flat.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "metazoa1k_rd_sparse" \
     -d afterok:$(get_jobid metazoa1k_rd_2) \
     -o ${DIR}/logs/rd_sparse.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "metazoa1k_rd_disk" \
     -o ${DIR}/logs/rd_disk.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 600 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "metazoa1k_rd_brwt" \
     -d afterok:$(get_jobid metazoa1k_rd_2) \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 5-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 10"

sbatch -J "metazoa1k_rd_brwt_relax" \
     -d afterok:$(get_jobid metazoa1k_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=18G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"
```



```bash
######################## Metazoa chunks ##########################

for i in {0,9}; do
    DIR=~/metagenome/data/cloudcompute/metazoa_graphs/nobackup/chunk_${i};
    METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
    mkdir -p ${DIR}/logs;

    cut -f1 ${DIR}/metazoa_all_cleanable_no_pacbio_no_nanopore_metadata_only_genomic_${i}_chunk.tsv | tail -n +2 > $DIR/ids.txt;
    for x in $(cat $DIR/ids.txt); do echo /cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa/clean/${x:0:3}/${x:0:6}/${x}/$x.fasta.gz; done > $DIR/samples.txt;

    sbatch -J "metazoa_${i}" \
         -o $DIR/logs/build_graph.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=19G \
         --partition=compute \
         --exclude compute-biomed-10 \
        --wrap="cat ${DIR}/samples.txt \
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
                $DIR/graph.dbg \
            && rm $DIR/graph.dbg; \
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

    mkdir ${DIR}/columns;

    sbatch -J "metazoa_${i}_annotate" \
         -d afterok:$(get_jobid metazoa_${i}) \
         -o $DIR/logs/annotate_graph.slog \
         -t 7-00 \
         --cpus-per-task 56 \
         --mem-per-cpu=15G \
         --exclude compute-biomed-10 \
        --wrap="cat ${DIR}/samples.txt \
            | /usr/bin/time -v $METAGRAPH annotate -v \
                -i $DIR/graph_primary.dbg \
                --anno-filename \
                --separately \
                -o ${DIR}/columns \
                -p 4 \
                --threads-each 14"

    rm -rf ${DIR}/rd;
    mkdir ${DIR}/rd;

    ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

    mkdir ${DIR}/rd/rd_columns;

    sbatch -J "metazoa_${i}_rd_0" \
         -d afterok:$(get_jobid metazoa_${i}_annotate) \
         -o $DIR/logs/rd_0.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=19G \
         --exclude compute-biomed-10 \
        --wrap="find ${DIR}/columns -name \"*.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 0 \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/rd/rd_columns/out \
                --mem-cap-gb 500 \
                -p 34"


    sbatch -J "metazoa_${i}_rd_1" \
         -d afterok:$(get_jobid metazoa_${i}_rd_0) \
         -o $DIR/logs/rd_1.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=19G \
         --exclude compute-biomed-10 \
        --wrap="find ${DIR}/columns -name \"*.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 1 \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/rd/rd_columns/out \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                --mem-cap-gb 500 \
                -p 34 && \
            rm ${DIR}/rd/rd_columns/*.row_count"

    sbatch -J "metazoa_${i}_rd_2" \
         -d afterok:$(get_jobid metazoa_${i}_rd_1) \
         -o $DIR/logs/rd_2.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=19G \
         --exclude compute-biomed-10 \
        --wrap="find ${DIR}/columns -name \"*.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --row-diff-stage 2 \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/rd/rd_columns/out \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                --mem-cap-gb 500 \
                -p 34 && \
            rm ${DIR}/rd/rd_columns/*.row_reduction && \
            rm ${DIR}/rd/graph.dbg.pred* && \
            rm ${DIR}/rd/graph.dbg.succ* && \
            ls -l ${DIR}/columns | grep annodbg > ${DIR}/columns.txt && \
            rm -r ${DIR}/columns"


    sbatch -J "metazoa_${i}_rd_brwt" \
         -d afterok:$(get_jobid metazoa_${i}_rd_2) \
         -o ${DIR}/logs/rd_brwt.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=24G \
        --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_brwt \
                --greedy \
                -i ${DIR}/rd/graph.dbg \
                -o ${DIR}/annotation \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -p 34 --parallel-nodes 4"

    sbatch -J "metazoa_${i}_rd_brwt_relax" \
         -d afterok:$(get_jobid metazoa_${i}_rd_brwt) \
         -o ${DIR}/logs/rd_brwt_relax.slog \
         -t 2-00 \
         --cpus-per-task 17 \
         --mem-per-cpu=40G \
        --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 17 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt.annodbg && \
        rm ${DIR}/annotation.row_diff_brwt.annodbg"
done


for i in 00 {0..15}; do
    DIR=~/metagenome/data/cloudcompute/metazoa_graphs/nobackup/chunk_${i};
    METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
    mkdir -p ${DIR}/logs;

    sbatch -J "stats_${i}" \
         -o $DIR/logs/graph_stats.slog \
         -t 7-00 \
         --cpus-per-task 34 \
         --mem-per-cpu=10G \
         --partition=compute \
        --wrap="/usr/bin/time -v $METAGRAPH stats -v --count-dummy \
                -p 34 \
                $DIR/graph_primary_small.dbg"
done
```



```bash

compute-biomed-16:~$ /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph query --query-mode matches         --num-top-labels 10 --min-kmers-fraction-label 0 -v         -i /cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/graph_primary.dbg         -a /cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/annotation.row_diff_flat.annodbg         /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/100_studies.fq > /dev/null

[2023-06-02 16:48:12.472] [trace] Query graph constructed for batch of sequences with 3302993 bases from '/cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/100_studies.fq' in 24.73789 sec, query redundancy: 1.38 bp/kmer, queried in 3.11557 sec
[2023-06-02 16:48:12.584] [trace] File '/cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/100_studies.fq' with 3302993 base pairs was processed in 27.965113572 sec, throughput: 118111.2 bp/s



cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=100_studies.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/graph_primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/annotation.row_diff_flat.annodbg;

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=microbe.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/BIGSI/graph_primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/BIGSI/annotation/annotation_primary_renamed.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=fungi.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/graph_primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/annotation_new.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=mouse.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/mus_musculus/metazoa_mus_musculus_primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/mus_musculus/nobackup/build/annotation.relaxed.relabeled.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=metazoa_1k.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_1k_studies/graph_primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_1k_studies/annotation.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=metasub.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/metasub/graphs/output_k41_cleaned_graph/graph_merged_k41.primary.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/metasub/graphs/k41/build/annotation.relaxed.row_diff_brwt.annodbg

sbatch -J "${FILE}.query" \
     -o ${FILE}.query.slog \
     -t 2-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH query -v --query-mode matches \
            --num-top-labels 10 --min-kmers-fraction-label 0 \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"

sbatch -J "${FILE}.query_align" \
     -o ${FILE}.query_align.slog \
     -t 2-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH query -v --query-mode matches \
            --align --num-top-labels 10 --min-kmers-fraction-label 0 \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"

sbatch -J "${FILE}.align" \
     -o ${FILE}.align.slog \
     -t 2-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH align -v \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"



METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
sbatch -J "convert" \
     -o /dev/null \
     -t 24-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=300G \
    --wrap="$METAGRAPH transform --state stat -o /cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/graph_primary.dbg /cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/graph_primary.small.dbg -p 1 -v"





cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=100_studies.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/graph_primary_small.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/annotation.row_diff_flat.annodbg;

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=microbe.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/BIGSI/graph_primary.small.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/BIGSI/annotation/annotation_primary_renamed.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=fungi.fq;
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/graph_primary.small.dbg;
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/fungi_graphs/annotation_new.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=mouse.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/mus_musculus/metazoa_mus_musculus_primary.small.dbg
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/mus_musculus/nobackup/build/annotation.relaxed.relabeled.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=metazoa_1k.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_1k_studies/graph_primary.small.dbg
ANNO=/cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa_graphs/nobackup/build_1k_studies/annotation.relaxed.row_diff_brwt.annodbg

cd ~/projects/projects2014-metagenome/metagraph/build_test;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
FILE=metasub.fq
GRAPH=/cluster/work/grlab/projects/metagenome/data/metasub/graphs/output_k41_cleaned_graph/graph_merged_k41.primary.small.dbg
ANNO=/cluster/work/grlab/projects/metagenome/data/metasub/graphs/k41/build/annotation.relaxed.row_diff_brwt.annodbg

sbatch -J "small.${FILE}.query" \
     -o small.${FILE}.query.slog \
     -t 24-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH query -v --query-mode matches \
            --num-top-labels 10 --min-kmers-fraction-label 0 \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"

sbatch -J "small.${FILE}.query_align" \
     -o small.${FILE}.query_align.slog \
     -t 24-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH query -v --query-mode matches \
            --align --num-top-labels 10 --min-kmers-fraction-label 0 \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"

sbatch -J "small.${FILE}.align" \
     -o small.${FILE}.align.slog \
     -t 24-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=400G \
    --wrap="/usr/bin/time -v $METAGRAPH align -v \
            -i $GRAPH -a $ANNO \
            /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/$FILE > /dev/null"




sbatch -J "query" \
     -o query.slog \
     -t 2-00 \
     --cpus-per-task 1 \
     --mem-per-cpu=60G \
    --wrap="/usr/bin/time -v $METAGRAPH query --query-mode matches \
        --num-top-labels 10 --min-kmers-fraction-label 0 -v \
        -i /cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/graph_primary.dbg \
        -a /cluster/work/grlab/projects/metagenome/data/cloudcompute/random_100_studies/annotation.row_diff_flat.annodbg \
        /cluster/work/grlab/projects/metagenome/data/cloudcompute/benchmark/queries/100_studies.fq"
```





```bash
######################## Random (100 studies) ##########################

DIR=~/metagenome/data/cloudcompute/random_100_studies;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

find ~/metagenome/data/cloudcompute/metagraph_1kstudies -name "*.fasta.gz" > $DIR/samples.txt;

mkdir -p ${DIR}/logs;

sbatch -J "random100" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=3G \
    --wrap="cat ${DIR}/samples.txt \
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
            $DIR/graph.dbg; \
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


mkdir ${DIR}/columns;

sbatch -J "random100_annotate" \
     -d afterok:$(get_jobid random100) \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=3G \
    --wrap="cat ${DIR}/samples.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 5 \
            --threads-each 8"


rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "random100_rd_0" \
     -d afterok:$(get_jobid random100_annotate) \
     -o $DIR/logs/rd_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 70 \
            -p 34"


sbatch -J "random100_rd_1" \
     -d afterok:$(get_jobid random100_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 70 \
            -p 34"


sbatch -J "random100_rd_2" \
     -d afterok:$(get_jobid random100_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 70 \
            -p 34"


sbatch -J "random100_rd_flat" \
     -d afterok:$(get_jobid random100_rd_2) \
     -o ${DIR}/logs/rd_flat.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random100_rd_sparse" \
     -d afterok:$(get_jobid random100_rd_2) \
     -o ${DIR}/logs/rd_sparse.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random100_rd_brwt" \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --subsample 10000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 10"

sbatch -J "random100_rd_brwt_relax" \
     -d afterok:$(get_jobid random100_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"

sbatch -J "random100_rd_disk" \
     -o ${DIR}/logs/rd_disk.slog \
     -t 00-12 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_disk \
            --mem-cap-gb 250 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"




sbatch -J "random100_rd_brwt_orig_linkage" \
     -o ${DIR}/logs/rd_brwt_orig_linkage.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="for x in \$(cat samples.txt); do echo columns/\$(basename \$x).column.annodbg; done \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --linkage \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/linkage_from_columns \
            -p 34"

sbatch -J "random100_rd_brwt_orig" \
     -d afterok:$(get_jobid random100_rd_brwt_orig_linkage) \
     -o ${DIR}/logs/rd_brwt_orig.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="for x in \$(cat samples.txt); do echo rd/rd_columns/\$(basename \$x).row_diff.annodbg; done \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --linkage-file ${DIR}/linkage_from_columns \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation_orig \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 10"


sbatch -J "random100_rd_brwt_orig_relax" \
     -d afterok:$(get_jobid random100_rd_brwt_orig) \
     -o ${DIR}/logs/rd_brwt_orig_relax_orig.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation_orig.relaxed \
        ${DIR}/annotation_orig.row_diff_brwt.annodbg"

```





```bash
DIR=~/metagenome/data/cloudcompute/random10K_cleaned_index/data;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir -p $DIR/logs;

find ~/metagenome/data/cloudcompute/random10K_unfiltered/ -name "*.fasta.gz" > $DIR/list.txt;

sbatch -J "random10K_cleaned_2" \
     -o $DIR/logs/clean_3.%A_%a.slog \
     -t 00-120 \
     --array=1-$(cat $DIR/list_3.txt | wc -l) \
     --cpus-per-task 4 \
     --mem-per-cpu=20G \
    --wrap="bash -c \"echo START \\\${SLURM_ARRAY_TASK_ID}; id=\\\$(sed -n \\\${SLURM_ARRAY_TASK_ID}p $DIR/list_3.txt); \
        echo \\\${id}; \
        file=/cluster/home/mikhaika/metagenome/data/cloudcompute/random10K_unfiltered/clean/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/\\\${id}.fasta.gz; \
        echo $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}; \
        rm -rf $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}; \
        mkdir -p $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}; \
        /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            --mode canonical \
            --inplace \
            --mem-cap-gb 20 \
            --count-kmers \
            --disk-swap ~/metagenome/scratch/nobackup \
            -p 4 \
            -o $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph \
            \\\${file} \
            > $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/build.log 2>&1; \
    /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --prune-tips 62 --prune-unitigs 0 --fallback 3 \
            -p 4 \
            -o $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/contigs_fb3 \
            $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph.dbg \
            > $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/clean_fb3.log 2>&1; \
    /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --prune-tips 62 --prune-unitigs 0 --fallback 2 \
            -p 4 \
            -o $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/contigs_fb2 \
            $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph.dbg \
            > $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/clean_fb2.log 2>&1; \
    /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --prune-unitigs 0 --fallback 3 \
            -p 4 \
            -o $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/contigs_fb3_with_tips \
            $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph.dbg \
            > $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/clean_fb3_with_tips.log 2>&1; \
    /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --prune-unitigs 0 --fallback 2 \
            -p 4 \
            -o $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/contigs_fb2_with_tips \
            $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph.dbg \
            > $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/clean_fb2_with_tips.log 2>&1; \
    rm $DIR/\\\${id:0:3}/\\\${id:0:6}/\\\${id}/graph.dbg*;\""
```


```bash
######################## (5K RANDOM SUBSET CLEANED) ##########################

DIR=~/metagenome/data/cloudcompute/random5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir -p ${DIR}/logs

sbatch -J "random5K" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
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
            $DIR/graph.dbg; \
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


DIR=~/metagenome/data/cloudcompute/random5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir ${DIR}/columns;

sbatch -J "random5K_annotate" \
     -d afterok:$(get_jobid random5K) \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 8 \
            --threads-each 8"


rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "random5K_rd_0" \
     -d afterok:$(get_jobid random5K_annotate) \
     -o $DIR/logs/rd_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 150 \
            -p 34"


sbatch -J "random5K_rd_1" \
     -d afterok:$(get_jobid random5K_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 150 \
            -p 34"


sbatch -J "random5K_rd_2" \
     -d afterok:$(get_jobid random5K_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 150 \
            -p 34"

DIR=~/metagenome/data/cloudcompute/random5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

sbatch -J "random5K_rd_flat" \
     -d afterok:$(get_jobid random5K_rd_2) \
     -o ${DIR}/logs/rd_flat.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random5K_rd_sparse" \
     -d afterok:$(get_jobid random5K_rd_2) \
     -o ${DIR}/logs/rd_sparse.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random5K_rd_brwt" \
     -d afterok:$(get_jobid random5K_rd_2) \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 10"

sbatch -J "random5K_rd_brwt_relax" \
     -d afterok:$(get_jobid random5K_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=8G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"




######### with counts ##########

DIR=~/metagenome/data/cloudcompute/random5k_cleaned_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir ${DIR}/columns_counts;

sbatch -J "random5K_anno_counts" \
     -o $DIR/logs/annotate_graph_counts.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=10G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            --count-kmers \
            --count-width 16 \
            --mem-cap-gb 20 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -o ${DIR}/columns_counts \
            -p 5 \
            --threads-each 8"

rm -rf ${DIR}/rd_counts;
mkdir ${DIR}/rd_counts;

ln -s $DIR/graph_primary.dbg $DIR/rd_counts/graph.dbg;

mkdir ${DIR}/rd_counts/rd_columns

sbatch -J "random5K_rd_counts_0" \
     -d afterok:$(get_jobid random5K_anno_counts) \
     -o $DIR/logs/rd_counts_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_counts_1" \
     -d afterok:$(get_jobid random5K_rd_counts_0) \
     -o $DIR/logs/rd_counts_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_counts_2" \
     -d afterok:$(get_jobid random5K_rd_counts_1) \
     -o $DIR/logs/rd_counts_2.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 450 \
            -p 34"

sbatch -J "random5K_rd_counts_brwt" \
     -d afterok:$(get_jobid random5K_rd_counts_2) \
     -o ${DIR}/logs/rd_counts_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/rd_counts/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 4"

sbatch -J "random5K_rd_counts_brwt_relax" \
     -d afterok:$(get_jobid random5K_rd_counts_brwt) \
     -o ${DIR}/logs/rd_counts_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_int_brwt.annodbg"


```





```bash
######################## (5K RANDOM SUBSET UNCLEANED) ##########################

DIR=~/metagenome/data/cloudcompute/random5K_unfiltered_singletons_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

mkdir -p ${DIR}/logs

sbatch -J "random5K" \
     -o $DIR/logs/build_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=25G \
    --wrap="cat ${DIR}/5k_list.txt \
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
            $DIR/graph.dbg; \
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


mkdir ${DIR}/columns_counts;

sbatch -J "random5K_anno_counts" \
     -d afterok:$(get_jobid random5K) \
     -o $DIR/logs/annotate_graph_counts.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            --count-kmers \
            --count-width 16 \
            --mem-cap-gb 20 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -o ${DIR}/columns_counts \
            -p 5 \
            --threads-each 8"

mkdir ${DIR}/columns;

sbatch -J "random5K_anno" \
     -d afterok:$(get_jobid random5K) \
     -o $DIR/logs/annotate_graph.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="cat ${DIR}/5k_list.txt \
        | /usr/bin/time -v $METAGRAPH annotate -v \
            -i $DIR/graph_primary.dbg \
            --anno-filename \
            --separately \
            -o ${DIR}/columns \
            -p 5 \
            --threads-each 8"


rm -rf ${DIR}/rd;
mkdir ${DIR}/rd;

ln -s $DIR/graph_primary.dbg $DIR/rd/graph.dbg;

mkdir ${DIR}/rd/rd_columns

sbatch -J "random5K_rd_0" \
     -d afterok:$(get_jobid random5K_anno) \
     -o $DIR/logs/rd_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --count-kmers \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_1" \
     -d afterok:$(get_jobid random5K_rd_0) \
     -o $DIR/logs/rd_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_2" \
     -d afterok:$(get_jobid random5K_rd_1) \
     -o $DIR/logs/rd_2.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 450 \
            -p 34"

sbatch -J "random5K_rd_flat" \
     -d afterok:$(get_jobid random5K_rd_2) \
     -o ${DIR}/logs/rd_flat.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_flat \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random5K_rd_sparse" \
     -d afterok:$(get_jobid random5K_rd_2) \
     -o ${DIR}/logs/rd_sparse.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_sparse \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 34"

sbatch -J "random5K_rd_brwt" \
     -o ${DIR}/logs/rd_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/rd/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 5"

sbatch -J "random5K_rd_brwt_relax" \
     -d afterok:$(get_jobid random5K_rd_brwt) \
     -o ${DIR}/logs/rd_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_brwt.annodbg"


######### with counts ##########

DIR=~/metagenome/data/cloudcompute/random5K_unfiltered_singletons_index;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;

rm -rf ${DIR}/rd_counts;
mkdir ${DIR}/rd_counts;

ln -s $DIR/graph_primary.dbg $DIR/rd_counts/graph.dbg;

mkdir ${DIR}/rd_counts/rd_columns

sbatch -J "random5K_rd_counts_0" \
     -d afterok:$(get_jobid random5K_anno) \
     -o $DIR/logs/rd_counts_0.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_counts_1" \
     -d afterok:$(get_jobid random5K_rd_counts_0) \
     -o $DIR/logs/rd_counts_1.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 650 \
            -p 34"


sbatch -J "random5K_rd_counts_2" \
     -d afterok:$(get_jobid random5K_rd_counts_1) \
     -o $DIR/logs/rd_counts_2.slog \
     -t 00-120 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="find ${DIR}/columns_counts -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --count-kmers \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/rd_counts/rd_columns/out \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            --mem-cap-gb 450 \
            -p 34"

sbatch -J "random5K_rd_counts_brwt" \
     -d afterok:$(get_jobid random5K_rd_counts_2) \
     -o ${DIR}/logs/rd_counts_brwt.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=24G \
    --wrap="find ${DIR}/rd_counts/rd_columns/ -name \"*.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy \
            -i ${DIR}/rd_counts/graph.dbg \
            -o ${DIR}/annotation \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 34 --parallel-nodes 4"

sbatch -J "random5K_rd_counts_brwt_relax" \
     -d afterok:$(get_jobid random5K_rd_counts_brwt) \
     -o ${DIR}/logs/rd_counts_brwt_relax.slog \
     -t 00-24 \
     --cpus-per-task 34 \
     --mem-per-cpu=19G \
    --wrap="/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -p 34 \
        --relax-arity 32 \
        -o ${DIR}/annotation.relaxed \
        ${DIR}/annotation.row_diff_int_brwt.annodbg"



```bash
sbatch -J "compute_metazoa_stats_SRR" \
     -o compute_stats_SRR.slog \
     -t 05-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=1G \
    --wrap="find clean/SRR -name \"*.fasta.gz\" | xargs -P 34 -n 1 ./compute_stats.sh > stats_SRR.txt"


cd /cluster/work/grlab/projects/metagenome/raw_data/refseq/release97/fna;
sbatch -J "compute_stats" \
     -o compute_stats.slog \
     -t 05-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=1G \
    --wrap="find complete_split -name \"*.fna.gz\" | xargs -P 34 -n 1 /cluster/work/grlab/projects/metagenome/data/cloudcompute/metazoa/compute_stats.sh > stats_refseq.txt"
```


######################## TCGA with counts ##########################

```bash
DIR=~/metagenome/data/tcga_counts/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
KMC=~/projects/projects2014-metagenome/metagraph/build_release/KMC/kmc;

mkdir $DIR;
mkdir $DIR/logs;
mkdir $DIR/graphs;
mkdir $DIR/unitigs;

find ~/metagenome/raw_data/tcga/data/ -name "*.r1.fq.gz" > $DIR/input_fastq.txt

sbatch -J "build_graphs" \
     --array=1-11095 \
     -o $DIR/logs/build_graphs.slog \
     -t 2-00 \
     --cpus-per-task 8 \
     --mem-per-cpu=15G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="file=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${DIR}/input_fastq.txt); \
        file=\${file%.r1.fq.gz}; \
        name=\$(basename \${file}); \
        [ -f $DIR/unitigs/\$name.fasta.gz ] && exit 0; \
        echo \"\${file}.r1.fq.gz\" > $DIR/graphs/\${name}.input; \
        [ -f \${file}.r2.fq.gz ] && echo \"\${file}.r2.fq.gz\" >> $DIR/graphs/\${name}.input; \
        mkdir ~/metagenome/scratch/nobackup/stripe_1/\${name}.kmc_cache; \
        /usr/bin/time -v $KMC -k31 -m40 -sm -ci1 -cs65535 -fq -t8 \
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
            -p 8 \
            -o $DIR/graphs/\$name \
            $DIR/graphs/\$name.kmc_suf > $DIR/graphs/\$name.build.log 2>&1 \
        && rm $DIR/graphs/\$name.kmc_* \
        && /usr/bin/time -v $METAGRAPH clean \
            --to-fasta --primary-kmers \
            --prune-unitigs 0 \
            --fallback 2 \
            -p 8 \
            -o $DIR/unitigs/\$name \
            $DIR/graphs/\$name.dbg > $DIR/graphs/\$name.clean.log 2>&1 \
        && rm $DIR/graphs/\$name.dbg*"


DIR=~/metagenome/data/tcga_counts/build;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;

sbatch -J "build_joint_graph" \
     -d afterok:2535396 \
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
DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
mkdir ${DIR};
mkdir ${DIR}/logs;
mkdir ${DIR}/graphs;
mkdir ${DIR}/unitigs;

sbatch -J "build_clean_graphs" \
     --array=1-11095 \
     -d afterok:2535396 \
     -o $DIR/logs/build_clean_graphs.slog \
     -t 2-00 \
     --cpus-per-task 4 \
     --mem-per-cpu=20G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="file=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${DIR}/../input_fastq.txt); \
        file=\${file%.r1.fq.gz}; \
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


WINDOW_SIZE=5;
DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
cd $DIR;
mkdir -p batches;
cd batches;
split -d -n r/40 <(find $DIR/unitigs -name "*.fasta.gz" | shuf);
mkdir -p ${DIR}/columns;

DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
for N in {0..39}; do
    N=$(printf "%02d" $N);
    list=x$N;
    sbatch -J "tcga_annotate_${WINDOW_SIZE}_${list}" \
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


DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/../graph_primary.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
sbatch -J "tcga_count_${WINDOW_SIZE}_rd" \
     -o ${DIR}/logs/count_rd.slog \
     -t 7-00 \
     --cpus-per-task 34 \
     --mem-per-cpu=13G \
     --partition=compute,gpu \
     --exclude gpu-biomed-12,gpu-biomed-23,gpu-biomed-10 \
    --wrap="find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb 400 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34 && \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb 400 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34  && \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb 400 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 34";

DIR=~/metagenome/data/tcga_counts/build/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
sbatch -J "tcga_count_${WINDOW_SIZE}_rd_brwt" \
     -d afterok:$(get_jobid tcga_count_${WINDOW_SIZE}_rd) \
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



######################## QUERY ##########################

```bash
~/time -v ./metagraph_DNA query -v --query-mode matches --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> query.fast && \
~/time -v ./metagraph_DNA query -v --query-mode matches --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> query.small_flat && \
~/time -v ./metagraph_DNA query -v --query-mode matches --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_brwt.annodbg -v sra/100_studies.fq > /dev/null 2> query.small && \
~/time -v ./metagraph_DNA query -v --query-mode matches --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> query.fast_disk && \
~/time -v ./metagraph_DNA query -v --query-mode matches --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> query.small_disk && \


~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> query_align_new.fast &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> query_align_new.small_flat &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_brwt.annodbg -v sra/100_studies.fq > /dev/null 2> query_align_new.small &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> query_align_new.fast_disk &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> query_align_new.small_disk &


~/time -v ./metagraph_DNA align -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> align_new.fast && \
~/time -v ./metagraph_DNA align -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v sra/100_studies.fq > /dev/null 2> align_new.small_flat && \
~/time -v ./metagraph_DNA align -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_brwt.annodbg -v sra/100_studies.fq > /dev/null 2> align_new.small && \
~/time -v ./metagraph_DNA align -i /optane/random_100_studies.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> align_new.fast_disk && \
~/time -v ./metagraph_DNA align -i /optane/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v sra/100_studies.fq > /dev/null 2> align_new.small_disk




~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i sra/random_100_studies.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v <(head -n 7400 sra/100_studies.fq) > /dev/null 2> query_align.fast &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i sra/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_flat.annodbg -v <(head -n 7400 sra/100_studies.fq) > /dev/null 2> query_align.small_flat &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i sra/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_brwt.annodbg -v <(head -n 7400 sra/100_studies.fq) > /dev/null 2> query_align.small &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i sra/random_100_studies.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v <(head -n 7400 sra/100_studies.fq) > /dev/null 2> query_align.fast_disk &
~/time -v ./metagraph_DNA query -v --query-mode matches --align --num-top-labels 10 --min-kmers-fraction-label 0 -i sra/random_100_studies.small.dbg -a /optane/random_100_studies.row_diff_disk.annodbg -v <(head -n 7400 sra/100_studies.fq) > /dev/null 2> query_align.small_disk &
```
