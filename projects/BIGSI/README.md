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
                --reverse \
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
