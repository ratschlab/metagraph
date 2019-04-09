# Graph Scalability Experiments

```bash
k=17; bsub -J "count_kmers[1-2789]%500" \
           -W 72:00 -n 10 -R "rusage[mem=3400]" \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/count_kmers_k$k.lsf
    "./count_kmers.sh $k \
        \"\$(sed -n \${LSB_JOBINDEX}p refseq.txt)\" \
        ~/metagenome/data/refseq/refseq_var_k/refseq_k$k"
```

### Build graphs and extract sequences

```bash
k=17; bsub -J contigs_k$k[1-2789]%100 \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/kmc_to_contigs_$k.lsf \
           -W 20:00 -n 12 -R "rusage[mem=10000] span[hosts=1]" \
    "x=\"\$(sed -n \${LSB_JOBINDEX}p refseq.txt | xargs -n 1 basename)\"; \
        ../../metagraph/scripts/kmc_to_contigs.sh \
            ~/metagenome/data/refseq/refseq_var_k/refseq_k$k/\$x.kmc_suf $k"

k=29; bsub -J contigs_k$k[1-2789]%100 \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/kmc_to_contigs_$k.lsf \
           -W 20:00 -n 12 -R "rusage[mem=10000] span[hosts=1]" \
    "x=\"\$(sed -n \${LSB_JOBINDEX}p refseq.txt | xargs -n 1 basename)\"; \
        ../../metagraph/scripts/kmc_to_contigs.sh \
            ~/metagenome/data/refseq/refseq_var_k/refseq_k$k/\$x.kmc_suf $k"

k=35; bsub -J contigs_k$k[1-2789]%100 \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/kmc_to_contigs_$k.lsf \
           -W 20:00 -n 12 -R "rusage[mem=12000] span[hosts=1]" \
    "x=\"\$(sed -n \${LSB_JOBINDEX}p refseq.txt | xargs -n 1 basename)\"; \
        ../../metagraph/scripts/kmc_to_contigs.sh \
            ~/metagenome/data/refseq/refseq_var_k/refseq_k$k/\$x.kmc_suf $k"

k=71; bsub -J contigs_k$k[1-599]%100 \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/kmc_to_contigs_$k.lsf \
           -W 20:00 -n 18 -R "rusage[mem=15000] span[hosts=1]" \
    "x=\"\$(sed -n \${LSB_JOBINDEX}p refseq.txt | xargs -n 1 basename)\"; \
        ../../metagraph/scripts/kmc_to_contigs.sh \
            ~/metagenome/data/refseq/refseq_var_k/refseq_k$k/\$x.kmc_suf $k"

k=71; bsub -J contigs_k$k[600-2789]%100 \
           -o $HOME/metagenome/data/refseq/refseq_var_k/refseq_k$k/kmc_to_contigs_$k.lsf \
           -W 20:00 -n 3 -R "rusage[mem=50000] span[hosts=1]" \
    "x=\"\$(sed -n \${LSB_JOBINDEX}p refseq.txt | xargs -n 1 basename)\"; \
        ../../metagraph/scripts/kmc_to_contigs.sh \
            ~/metagenome/data/refseq/refseq_var_k/refseq_k$k/\$x.kmc_suf $k"
```

### Build graphs for nested subsets

#### Succinct

```bash
k=17; for i in {1..90}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}.lsf \
         -W 20:00 -n 18 -R "rusage[mem=36700] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --mem-cap-gb 600 -k $k -p 36 \
            -o $data_root/subgraphs_k$k/subset_${i}"; \
done


k=29; for i in {1..55}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}.lsf \
         -W 20:00 -n 18 -R "rusage[mem=13100] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --mem-cap-gb 220 -k $k -p 36 -s 40 \
            -o $data_root/subgraphs_k$k/subset_${i}"; \
done


k=35; for i in {1..55}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}.lsf \
         -W 20:00 -n 18 -R "rusage[mem=13100] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --mem-cap-gb 220 -k $k -p 36 -s 40 \
            -o $data_root/subgraphs_k$k/subset_${i}"; \
done


k=71; for i in {1..50}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}.lsf \
         -W 20:00 -n 18 -R "rusage[mem=16500] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --mem-cap-gb 280 -k $k -p 36 -s 40 \
            -o $data_root/subgraphs_k$k/subset_${i}"; \
done
```

#### Bitmap

```bash
k=17; for i in {1..90}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_bitmap_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_bitmap.lsf \
         -W 20:00 -n 18 -R "rusage[mem=13000] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph bitmap \
            --mem-cap-gb 180 -k $k -p 36 -s 4 \
            -o $data_root/subgraphs_k$k/subset_${i}_bitmap"; \
done


k=29; for i in {1..55}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_bitmap_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_bitmap.lsf \
         -W 20:00 -n 18 -R "rusage[mem=19400] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph bitmap \
            --mem-cap-gb 180 -k $k -p 36 -s 16 \
            -o $data_root/subgraphs_k$k/subset_${i}_bitmap"; \
done
```

#### Hash-based

```bash
k=17; for i in {1..20}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_hash_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_hash.lsf \
         -W 20:00 -n 1 -R "rusage[mem=$((25000 * i))] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph hash \
            -k $k \
            -o $data_root/subgraphs_k$k/subset_${i}_hash"; \
done


k=29; for i in {1..20}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_hash_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_hash.lsf \
         -W 20:00 -n 1 -R "rusage[mem=$((30000 * i))] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph hash \
            -k $k \
            -o $data_root/subgraphs_k$k/subset_${i}_hash"; \
done


k=35; for i in {1..20}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_hash_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_hash.lsf \
         -W 20:00 -n 1 -R "rusage[mem=$((40000 * i))] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph hash \
            -k $k \
            -o $data_root/subgraphs_k$k/subset_${i}_hash"; \
done


k=71; for i in {1..20}; do \
    data_root="$HOME/metagenome/data/refseq/refseq_var_k"; \
    bsub -J build_hash_k${k}_${i} \
         -oo $data_root/subgraphs_k$k/subset_${i}_hash.lsf \
         -W 20:00 -n 1 -R "rusage[mem=$((50000 * i))] span[hosts=1]" \
        "for line in \$(cat refseq_subsets/subset_$i); do \
            echo $data_root/refseq_k$k/\${line}.contigs.fasta.gz; \
        done \
        | gtime -v ../../metagraph/build/metagengraph build -v \
            --graph hash \
            -k $k \
            -o $data_root/subgraphs_k$k/subset_${i}_hash"; \
done
```

## Measure graph stats

```bash
ls -1aS $HOME/metagenome/data/refseq/refseq_var_k/subgraphs_k*/subset_*.dbg \
    | xargs -n 1 -I % bash -c \
        'memcap="$(($(du --block-size=1M % | cut -f1) * 5 / 4 + 100))"; \
        bsub -J "stats_succinct" -W 3:00 -n 1 -R "rusage[mem=$memcap] span[hosts=1]" \
             -oo %.stats \
            "gtime -v ../../metagraph/build/metagengraph stats %"'
```

```bash
ls -1aS $HOME/metagenome/data/refseq/refseq_var_k/subgraphs_k*/subset_*.bitmapdbg \
    | xargs -n 1 -I % bash -c \
        'memcap="$(($(du --block-size=1M % | cut -f1) * 5 / 4 + 100))"; \
        bsub -J "stats_bitmap" -W 3:00 -n 1 -R "rusage[mem=$memcap] span[hosts=1]" \
             -oo %.stats \
            "gtime -v ../../metagraph/build/metagengraph stats %"'
```

```bash
ls -1aS $HOME/metagenome/data/refseq/refseq_var_k/subgraphs_k*/subset_*.orhashdbg \
    | xargs -n 1 -I % bash -c \
        'memcap="$(($(du --block-size=1M % | cut -f1) * 6))"; \
        bsub -J "stats_hash" -W 3:00 -n 1 -R "rusage[mem=$memcap] span[hosts=1]" \
             -oo %.stats \
            "gtime -v ../../metagraph/build/metagengraph stats %"'
```


## Dependence on k

### Extract k-mers

```bash
for k in {10..85}; do \
    outdir="$HOME/metagenome/data/refseq/refseq_var_k/complete.1161.1.genomic.fna.gz/complete.1161.1.genomic.fna.gz.k$k"; \
    mkdir $outdir; \
    bsub -J "count_kmers" \
         -oo $outdir/complete.1161.1.genomic.fna.gz.k$k.count_kmers.lsf \
         -W 72:00 -n 2 -R "rusage[mem=8000]" \
        "./count_kmers.sh $k \
            ~/metagenome/raw_data/refseq/fna/complete.1161.1.genomic.fna.gz \
            $outdir"; \
done
```

### Build graphs and extract contigs

```bash
for k in {10..42}; do \
    dir="$HOME/metagenome/data/refseq/refseq_var_k/complete.1161.1.genomic.fna.gz/complete.1161.1.genomic.fna.gz.k$k"; \
    bsub -J "contigs" \
         -oo $dir/complete.1161.1.genomic.fna.gz.k$k.kmc_to_contigs.lsf \
         -W 72:00 -n 1 -R "rusage[mem=85000]" \
        "../../metagraph/scripts/kmc_to_contigs.sh \
            $outdir/complete.1161.1.genomic.fna.gz.kmc_suf $k";\
done
```

```bash
for k in {43..85}; do \
    dir="$HOME/metagenome/data/refseq/refseq_var_k/complete.1161.1.genomic.fna.gz/complete.1161.1.genomic.fna.gz.k$k"; \
    bsub -J "contigs" \
         -oo $dir/complete.1161.1.genomic.fna.gz.k$k.kmc_to_contigs.lsf \
         -W 72:00 -n 5 -R "rusage[mem=50000]" \
        "../../metagraph/scripts/kmc_to_contigs.sh \
            $outdir/complete.1161.1.genomic.fna.gz.kmc_suf $k";\
done
```
