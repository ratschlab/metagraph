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
```

## Annotate
```bash
mkdir temp
find ~/metagenome/data/BIGSI/ -name "*fasta.gz" > temp/files_to_annotate.txt
cd temp; split -l 5000 files_to_annotate.txt; cd ..

cd temp
for list in x*; do
    bsub -J "annotate_${list}" \
         -oo ~/metagenome/data/BIGSI/annotate_${list}.lsf \
         -W 96:00 \
         -n 15 -R "rusage[mem=23000] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build/metagraph annotate -v \
                -i ~/metagenome/data/BIGSI/graph \
                --parallel 15 \
                --anno-filename \
                --separately \
                2>&1 | tee ~/metagenome/data/BIGSI/annotate_${list}.log"; \
done
cd ..
```

## Generate subsets
```bash
cat temp/files_to_annotate.txt | shuf > subsets/files_shuffled.txt
cd subsets
for i in {1..20}; do N=$((750 * i)); head -n $N files_shuffled.txt > files_$N.txt; done
```

### Build graphs for subsets
```bash
mkdir ~/metagenome/data/BIGSI/subsets
for i in {1..20}; do
    N=$((750 * i));
    bsub -J "subset_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.lsf \
         -W 24:00 \
         -n 15 -R "rusage[mem=23000] span[hosts=1]" \
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

## Annotate graphs for subsets
```bash
for i in {20..1}; do
    N=$((750 * i));
    mkdir ~/metagenome/data/BIGSI/subsets/graph_subset_${N}
    bsub -J "annotate_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/annotate_subset_${N}.lsf \
         -W 96:00 \
         -n 15 -R "rusage[mem=3000] span[hosts=1]" \
        "cat subsets/files_${N}.txt \
            | /usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA annotate -v \
                -i ~/metagenome/data/BIGSI/subsets/graph_subset_${N}.dbg \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/BIGSI/subsets/graph_subset_${N} \
                2>&1"; \
done

for i in {1..20}; do
    N=$((750 * i));
    bsub -J "merge_columns_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/merge_columns_subset_${N}.lsf \
         -W 96:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "find ~/metagenome/data/BIGSI/subsets/graph_subset_${N}/ -name *.annodbg | sort \
            | /usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA merge_anno -v \
                -o ~/metagenome/data/BIGSI/subsets/annotation_subset_${N} \
                --parallel 20 \
                2>&1"; \
done
```

## Transform annotation
```bash
for i in {1..20}; do
    N=$((750 * i));
    bsub -J "to_row_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/columns_to_rows_subset_${N}.lsf \
         -W 96:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "/usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA transform_anno -v \
                --anno-type row \
                -o ~/metagenome/data/BIGSI/subsets/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation_subset_${N}.column.annodbg \
                --parallel 20 \
                2>&1"; \
done


for i in {20..1}; do
    N=$((750 * i));
    bsub -J "to_flat_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/rows_to_flat_subset_${N}.lsf \
         -W 96:00 \
         -n 1 -R "rusage[mem=$((N * 10 + 15000))] span[hosts=1]" \
        "/usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA transform_anno -v \
                --anno-type flat \
                -o ~/metagenome/data/BIGSI/subsets/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation_subset_${N}.row.annodbg \
                2>&1"; \
done


for i in {1..20}; do
    N=$((750 * i));
    bsub -J "to_brwt_${N}" \
         -oo ~/metagenome/data/BIGSI/subsets/columns_to_brwt_subset_${N}.lsf \
         -W 96:00 \
         -n 10 -R "rusage[mem=${N}] span[hosts=1]" \
        "/usr/bin/time -v ~/metagenome/metagraph_server/metagraph_DNA transform_anno -v \
                --anno-type brwt --greedy \
                -o ~/metagenome/data/BIGSI/subsets/annotation_subset_${N} \
                ~/metagenome/data/BIGSI/subsets/annotation_subset_${N}.column.annodbg \
                --parallel 20 \
                2>&1"; \
done
```

## Generate BIGSI bloom filters
```bash
bsub -J "bigsi_bloom[1-15000]%300" \
     -o ~/metagenome/data/BIGSI/subsets/bigsi/build_bloom.lsf \
     -W 20:00 \
     -n 1 -R "rusage[mem=3000] span[hosts=1]" \
     "file=\"\$(sed -n \${LSB_JOBINDEX}p subsets/files_15000.txt)\"; \
        x=\"\$(echo \$file | xargs -n 1 basename)\"; \
        sra_id=\"\${x%.unitigs.fasta.gz}\"; \
        ctx_file=\"/cluster/work/grlab/projects/metagenome/raw_data/BIGSI/ctx_subsets/\${sra_id}.ctx\"; \
        /usr/bin/time -v bigsi bloom \
            --config=/cluster/home/mikhaika/stuff/BIGSI/berkleydb.yaml \
            \$ctx_file \
            /cluster/work/grlab/projects/metagenome/data/BIGSI/subsets/bigsi/bloom/\${sra_id}.bloom \
            2>&1"
```
