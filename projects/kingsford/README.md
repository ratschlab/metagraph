```bash
find ~/metagenome/data/kingsford/kmc_20/ -name "*kmc_suf" > ~/metagenome/data/kingsford/kingsford_list_kmc.txt

bsub -J "build_single[1-2652]%500" \
     -o ~/metagenome/data/kingsford/build_single.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
        "file=\"\$(sed -n \${LSB_JOBINDEX}p ~/metagenome/data/kingsford/kingsford_list_kmc.txt)\"; \
        /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA build -v \
            -k 20 \
            --canonical \
            --count-kmers \
            --mem-cap-gb 8 \
            -p 1 \
            -o ~/metagenome/data/kingsford/single_graphs/\$(basename \${file%.fasta.gz.kmc.kmc_suf}) \
            \$file \
        2>&1"

bsub -w "done(build_single[*])" \
     -J "extract_contigs[1-2652]%500" \
     -o ~/metagenome/data/kingsford/extract_contigs.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=2000] span[hosts=1]" \
        "file=\"\$(sed -n \${LSB_JOBINDEX}p ~/metagenome/data/kingsford/kingsford_list_kmc.txt)\"; \
        file=\$(basename \${file%.fasta.gz.kmc.kmc_suf})
        /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph_DNA clean -v \
            --to-fasta --primary-kmers \
            -p 1 \
            -o ~/metagenome/data/kingsford/single_graphs/\$file \
            ~/metagenome/data/kingsford/single_graphs/\${file}.dbg \
        2>&1"

bsub -J "kingsford_build" \
     -oo ~/metagenome/data/kingsford/build_canonical.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=3000] span[hosts=1]" \
    "find ~/metagenome/data/kingsford/single_graphs/ -name \"*.fasta.gz\" \
        | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
            -k 20 \
            --canonical \
            --mem-cap-gb 80 \
            -p 36 \
            -o ~/metagenome/data/kingsford/kingsford_canonical \
            2>&1"

bsub -J "kingsford_to_primary" \
     -oo ~/metagenome/data/kingsford/extract_primary_contigs.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=2000] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform -v \
            --to-fasta --primary-kmers \
            -o ~/metagenome/data/kingsford/kingsford_primary \
            ~/metagenome/data/kingsford/kingsford_canonical.dbg \
            -p 36 \
            2>&1"

bsub -J "kingsford_build" \
     -oo ~/metagenome/data/kingsford/build_primary.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1500] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
            -k 20 \
            --mem-cap-gb 40 \
            -p 36 \
            -o ~/metagenome/data/kingsford/kingsford_primary \
            ~/metagenome/data/kingsford/kingsford_primary.fasta.gz \
            2>&1"


mkdir temp
cd temp
find ~/metagenome/data/kingsford/single_graphs/ -name "*.fasta.gz" > files_to_annotate.txt
split -l 150 files_to_annotate.txt

for list in x*; do
    bsub -J "kingsford_annotate_${list}" \
         -oo ~/metagenome/data/kingsford/annotation/annotate_primary_${list}.lsf \
         -W 4:00 \
         -n 15 -R "rusage[mem=1000] span[hosts=1]" \
        "cat ${list} \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                -i ~/metagenome/data/kingsford/kingsford_primary.dbg \
                --canonical \
                --count-kmers \
                --parallel 15 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/kingsford/annotation/columns_primary \
                2>&1"; \
done
cd ..

find ~/metagenome/data/kingsford/annotation/columns_primary/ -name "*.column.annodbg" \
    > ~/metagenome/data/kingsford/annotation/columns_primary.txt

bsub -J "kingsford_build_rbbrwt_relax20_primary" \
     -oo ~/metagenome/data/kingsford/annotation/build_rbbrwt_primary.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=6300] span[hosts=1]" \
    "cat ~/metagenome/data/kingsford/annotation/columns_primary.txt \
        | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --anno-type rb_brwt --relax-arity 20 \
            -o ~/metagenome/data/kingsford/annotation/annotation_primary \
            -p 36 \
            2>&1"

bsub -J "kingsford_cluster_columns_primary" \
     -oo ~/metagenome/data/kingsford/annotation/cluster_columns_primary.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=2100] span[hosts=1]" \
    "cat ~/metagenome/data/kingsford/annotation/columns_primary.txt \
        | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --linkage --greedy \
            --subsample 100000000 \
            -o ~/metagenome/data/kingsford/annotation/linkage_kingsford_primary.csv \
            -p 36 \
            2>&1";

bsub -J "kingsford_build_brwt_primary" \
     -oo ~/metagenome/data/kingsford/annotation/build_brwt_primary.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1000] span[hosts=1]" \
    "cat ~/metagenome/data/kingsford/annotation/columns_primary.txt \
        | gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
            --anno-type brwt \
            -i ~/metagenome/data/kingsford/annotation/linkage_kingsford_primary.csv \
            -o ~/metagenome/data/kingsford/annotation/annotation_primary \
            -p 36 --parallel-nodes 10 \
            2>&1"

bsub -J "kingsford_relax_brwt_primary" \
     -oo ~/metagenome/data/kingsford/annotation/relax_brwt_primary.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=1000] span[hosts=1]" \
    "gtime -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
        -p 36 \
        --relax-arity 16 \
        -o ~/metagenome/data/kingsford/annotation/annotation_primary.relaxed \
        ~/metagenome/data/kingsford/annotation/annotation_primary.brwt.annodbg \
        2>&1"
```

