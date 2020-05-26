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
                --canonical \
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
    "find ~/metagenome/data/BIGSI/annotation/columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph transform_anno -v \
            --linkage \
            --subsample 5000000 \
            -o ~/metagenome/data/BIGSI/annotation/linkage_BIGSI.csv \
            --parallel 96 \
            2>&1";

bsub -J "cluster" \
     -oo ~/metagenome/data/BIGSI/logs/cluster_columns_1M.lsf \
     -W 120:00 \
     -n 48 -R "rusage[mem=37500] span[hosts=1]" \
    "find ~/metagenome/data/BIGSI/annotation/columns/ -name \"*.column.annodbg\" \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_release/metagraph transform_anno -v \
            --linkage \
            --subsample 1000000 \
            -o ~/metagenome/data/BIGSI/annotation/linkage_BIGSI_1M.csv \
            --parallel 96 \
            2>&1";
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
                --canonical \
                --parallel 30 \
                --mem-cap-gb 300 \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N} \
                2>&1"; \
done
```

### Transform graphs
```bash
for i in {1..33}; do
    N=$((750 * i));
    bsub -J "to_small_graph_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/to_small_graph_${N}.lsf \
         -W 1:00 \
         -n 1 -R "rusage[mem=50000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform -v \
                --state small \
                --index-ranges 10 \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.small.indexed.dbg \
                ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.small.dbg \
                2>&1"; \
done
```


## Annotate graphs for subsets
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

extension=".column.annodbg";
for i in {33..1}; do
    N=$((750 * i));
    annotation="~/metagenome/data/BIGSI/subsets/annotation_subset_${N}${extension}";
    old="${annotation%$extension}.path_labels${extension}";
    mv $annotation $old;
    bsub -J "rename_columns_${N}" \
         -oo /dev/null \
         -W 1:00 \
         -n 1 -R "rusage[mem=$((N * 6))] span[hosts=1]" \
        "~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno \
            --rename-cols ~/metagenome/data/BIGSI/subsets/rename_columns.txt \
            -o $annotation \
            $old"; \
done
```

## Transform annotation
```bash
for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_row_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/column_to_rowy_${N}.lsf \
         -W 24:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type row \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.column.annodbg \
                --parallel 20 \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_flat_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/row_to_flat_${N}.lsf \
         -W 24:00 \
         -n 1 -R "rusage[mem=$((N * 10 + 15000))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type flat \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.row.annodbg \
                2>&1"; \
done


for i in {33..1}; do
    N=$((750 * i));
    bsub -J "to_rbfish_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/lsf_logs/row_to_rbfish_${N}.lsf \
         -W 20:00 \
         -n 1 -R "rusage[mem=$((N * 17))] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA transform_anno -v \
                --anno-type rbfish \
                -o ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${N}.row.annodbg \
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
```

## Query graph
```bash
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA

# file to query
for QUERY in ~/metagenome/data/BIGSI/subsets/query/samples/haib18CEM5453_HMCMJCCXY_SL336225.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/nucleotide_fasta_protein_homolog_model.fasta \
                ~/metagenome/data/BIGSI/subsets/query/samples/DRR067889.fasta; do
    NAME=metagraph.small_indexed.brwt_relax
    # name of the output folder
    OUTDIR=~/metagenome/data/BIGSI/subsets/query_results/$(basename $QUERY)/${NAME}
    mkdir -p $OUTDIR

    for num_columns in $(seq 750 3000 24750); do
        run="$METAGRAPH query -v --discovery-fraction 0.0 --count-labels --fast \
                -i \${TMPDIR}/graph.dbg \
                -a \${TMPDIR}/graph.brwt.annodbg \
                $QUERY"

        bsub -J "${NAME}.${num_columns}" \
            -W 3:50 \
            -n 36 -R "rusage[mem=4000] span[hosts=1] select[model==XeonGold_6140]" \
            -oo ${OUTDIR}/${num_columns}.lsf \
            " \
                TMPDIR=/dev/shm/$NAME_$num_columns_\${LSB_JOBID};
                mkdir \${TMPDIR};
                cp ~/metagenome/data/BIGSI/subsets/graph_subset_${num_columns}.dbg \
                    \${TMPDIR}/graph.dbg; \
                cp ~/metagenome/data/BIGSI/subsets/annotation/annotation_subset_${num_columns}.relaxed.brwt.annodbg \
                    \${TMPDIR}/graph.brwt.annodbg; \
                /usr/bin/time -v $run > /dev/null 2> /dev/null; \
                for i in {1..10}; do \
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
