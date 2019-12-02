

```bash
DIR='/cluster/work/grlab/projects/metagenome/benchmark_human_metagenome/nobackup/human_gut_sra'
bsub -J "trim[1-20000]%500" -W 72:00 -n 1 -R "rusage[mem=1000]" \
    "../trimming.sh ${DIR}/\"\$(sed -n \${LSB_JOBINDEX}p files_sorted_desc_5.txt)\" 0.02"
```
