# MetaSUB

## Sample processing - version 1 (see analysis_v1)

### Count k-mers with KMC

Count all k-mers using KMC3. 
```bash
bash run_count_kmers.sh
```
This will invoke
```bash
count_kmers.sh
```

Create a list of all kmc dumps
```bash
cat -1S ~/big_graph/metasub_wasabi/*.kmc_suf > metasub_kmc.txt
```

### Build graphs from KMC
Generate a graph from each kmc file:
```bash
bash run_kmc_to_graphs_array.sh
```

### Extract contigs from single graphs

```bash
bsub -J contigs_k15[1-4173]%600 \
    -o seq_extraction.lsf \
    -W 20:00 -n 1 -R "rusage[mem=22000] span[hosts=1]" \
    "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p metasub_kmc.txt)\"; \
        x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
        $METAGRAPH assemble \
            -o ~/big_graph/metasub/metasub_k15/\${x%.kmc_suf}.k15.sequences \
            \${kmc_file%.kmc_suf}.k15.dbg"
```

### Build graph
```bash
$METAGRAPH build --mode canonical --complete --graph bitmap \
    -o graph_17_complete_canonical.bitmapdbg -k 17
```

```bash
$METAGRAPH build --mode canonical --complete --graph bitmap \
    -o graph_15_complete_canonical.bitmapdbg -k 15
```

### Build annotation

#### Construct columns
Generate an annotation columns for all contigs of a sample
```bash
bash run_annotate_contigs.sh
```

#### Transform to BRWT

```bash
bsub -J convert_metasub_to_brwt_pm \
    -o conversion_to_brwt_f2_k13.lsf \
    -W 150:00 -n 30 -R "rusage[mem=10000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH transform_anno -v -p 60 \
        -o ~/big_graph/metasub/metasub_f2_k13_brwt_pm \
        --anno-type brwt --greedy \
        ~/big_graph/metasub_wasabi_graph_anno.k13.f2.column.annodbg \
        2>&1"
```

```bash
bsub -J convert_metasub_to_brwt_pm \
    -oo conversion_to_brwt.lsf \
    -W 150:00 -n 30 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH transform_anno -v \
        -o ~/big_graph/metasub/metasub_f2_k15_brwt_pm \
        --anno-type brwt --greedy -p 60 \
        ~/big_graph/metasub_wasabi_graph_anno.column.annodbg \
        2>&1"
```

#### Optimize BRWT

```bash
bsub -J relax_metasub_f2_k13_brwt_40 \
    -oo relax_brwt_f2_k13_40.lsf \
    -W 40:00 -n 5 -R "rusage[mem=4000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -o metasub_f2_k13_brwt_pm_40 \
        --relax-arity 40 -p 10 \
        metasub_f2_k13_brwt_pm.brwt.annodbg \
        2>&1"
```

```bash
for arity in {5,7,10,20,30,40}; \
    do bsub -J relax_metasub_f2_k15_brwt_${arity} \
        -oo relax_brwt_f2_k15_${arity}.lsf \
        -W 80:00 -n 5 -R "rusage[mem=50000] span[hosts=1]" \
        "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -o metasub_f2_k15_brwt_pm_${arity} \
            --relax-arity ${arity} -p 10 \
            metasub_f2_k15_brwt_pm.brwt.annodbg \
            2>&1"; \
done
```

```bash
bsub -J relax_metasub_f2_k17_brwt_40 \
    -oo relax_brwt_f2_k17_40.lsf \
    -W 100:00 -n 5 -R "rusage[mem=136000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
        -o metasub_f2_k17_brwt_pm_40 \
        --relax-arity 40 -p 10 \
        ../metasub_wasabi_graph_anno.k17.f2.brwt_pm.brwt.annodbg \
        2>&1"
```

#### Rename columns

```bash
$METAGRAPH transform_anno --rename-cols rename_columns.txt \
    -o metasub_f2_k17_brwt_pm_final.brwt.annodbg \
    ../metasub_wasabi_graph_anno.k17.f2.brwt_pm_40.brwt.annodbg
```

## Sample processing - version 2 (see analysis_v2)

## Start jobs with MetaSUB database

```bash
bsub -J metasub_graph_server15 \
    -oo /dev/null \
    -W 240:00 -n 10 -R "rusage[mem=19000] span[hosts=1]" \
    "$METAGRAPH server_query \
        -i graph_15_complete_canonical.bitmapdbg \
        -a metasub_f2_k15_brwt_pm_40_final.brwt.annodbg \
        -p 20 --port 42623 2>&1 \
        | tee metasub_graph_server.log"
```

```bash
bsub -J metasub_graph_server17 \
    -oo /dev/null \
    -W 240:00 -n 10 -R "rusage[mem=61000] span[hosts=1]" \
    "$METAGRAPH server_query \
        -i graph_17_complete_canonical.bitmapdbg \
        -a metasub_f2_k17_brwt_pm_20_final.brwt.annodbg \
        -p 20 --port 42623 2>&1 \
        | tee metasub_graph_server_k17.log"
```
