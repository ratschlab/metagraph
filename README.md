# Metagenome Graph Project

## Install

### Prerequisites
- cmake 3.10 or higher
- GNU GCC with C++17 (gcc-8.0.1 or higher), LLVM Clang (clang-7 or higher), or AppleClang (clang-1100.0.33.8 or higher)
- bzip2
- HTSlib

#### Optional
- boost and jemalloc-4.0.0 or higher (to build with *folly* for efficient small vector support)
- Python 3 (for running integration tests)

For compiling with **AppleClang**, the prerequisites can be installed as easy as:
```
brew install libomp cmake make bzip2 htslib boost jemalloc
```

For **Ubuntu** (20.04 LTS or higher) or **Debian** (10 or higher)
```
sudo apt-get install cmake libbz2-dev libhts-dev libjemalloc-dev libboost-all-dev
```

For **CentOS** (8 or higher)
```
yum install cmake bzip2-devel htslib-devel jemalloc-devel boost-devel
```

All prerequisites can also be installed by users **without root** using [brew](https://brew.sh) or [linuxbrew](https://linuxbrew.sh).


### Compile
1. `git clone --recursive https://github.com/ratschlab/metagraph.git`
2. install *sdsl-lite* in `metagraph/external-libraries/` by running the following script from the repository root directory
```bash
git submodule sync
git submodule update --init --recursive

pushd metagraph/external-libraries/sdsl-lite
./install.sh $PWD
popd
```

3. make a **build** directory `mkdir -p metagraph/build && cd metagraph/build`
4. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))` (for alphabets other than DNA, see below)
5. (optional) run unit tests `./unit_tests`
6. (optional) run integration tests `./integration_tests`

### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile|GProfile]` -- build modes (`Release` by default)
- `-DBUILD_STATIC=[ON|OFF]` -- link statically (`OFF` by default)
- `-DLINK_OPT=[ON|OFF]` -- enable link time optimization (`OFF` by default)
- `-DBUILD_KMC=[ON|OFF]` -- compile the KMC executable (`ON` by default)
- `-DWITH_AVX=[ON|OFF]` -- compile with support for the avx instructions (`ON` by default, if available)
- `-DWITH_MSSE42=[ON|OFF]` -- compile with support for the msse4.2 instructions (`ON` by default, if available)
- `-DCMAKE_DBG_ALPHABET=[Protein|DNA|DNA5|DNA_CASE_SENSITIVE]` -- alphabet to use (`DNA` by default)

## Typical workflow
1. Build de Bruijn graph from Fasta files, FastQ files, or [KMC k-mer counters](https://github.com/refresh-bio/KMC/):\
`./metagraph build`
2. Annotate graph using the column compressed annotation:\
`./metagraph annotate`
3. Transform the built annotation to a different annotation scheme:\
`./metagraph transform_anno`
4. Query annotated graph\
`./metagraph query`

### Example
```
DATA="../tests/data/transcripts_1000.fa"

./metagraph build -k 12 -o transcripts_1000 $DATA

./metagraph annotate -i transcripts_1000.dbg --anno-filename -o transcripts_1000 $DATA

./metagraph query -i transcripts_1000.dbg -a transcripts_1000.column.annodbg $DATA

./metagraph stats -a transcripts_1000.column.annodbg transcripts_1000.dbg
```

### Graph cleaning
```
DATA="../tests/data/transcripts_1000.fa"
K=12

./metagraph build -k $K --count-kmers -o transcripts_1000 $DATA

./metagraph clean --prune-tips $((2*K)) --prune-unitigs 0 --fallback 2 --to-fasta -o transcripts_1000_clean_contigs transcripts_1000.dbg

zless transcripts_1000_clean_contigs.fasta.gz | tail
```

For real examples, see [scripts](./metagraph/scripts).

### Print usage
`./metagraph`

### Build graph

* #### Simple build
```bash
./metagraph build -v --parallel 30 -k 20 --mem-cap-gb 10 \
                        -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz \
2>&1 | tee <LOG_DIR>/log.txt
```

* #### Build with disk swap (use to limit the RAM usage)
```bash
./metagraph build -v --parallel 30 -k 20 --mem-cap-gb 10 --disk-swap <GRAPH_DIR> \
                        -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz \
2>&1 | tee <LOG_DIR>/log.txt
```

* #### Build from chunks (use only for very large graphs)
1) Build chunks
```bash
for F in {\$,A,C,G,T,N}{\$,A,C,G,T,N}; do \
    ./metagraph build -v --parallel 30 -k 20 --mem-cap-gb 100 \
                            -o <GRAPH_DIR>/graph --suffix $F \
                            <DATA_DIR>/*.fasta.gz \
    2>&1 | tee <LOG_DIR>/log_$F.txt; \
done
```
2) Concatenate chunks
```bash
./metagraph concatenate -l 2 -i <GRAPH_DIR>/graph -o <GRAPH_DIR>/graph
```

#### Build from k-mers filtered with KMC
```bash
CUTOFF=5
K=20
./KMC/kmc -ci$CUTOFF -t30 -k$K -m5 -fq -b <FILE>.fasta.gz <FILE>.kmc_$CUTOFF ./KMC
./metagraph build -v -p 30 -k $K --mem-cap-gb 10 --kmc -o graph <FILE>.kmc_$CUTOFF
```

#### Distributed build
1) Build chunks
```bash
for F in {\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}; do \
    bsub -J assemble$F -W 8:00 -n 30 -R "rusage[mem=15000]" \
        "ls -1a <DATA_DIR>/*.fasta.gz | /usr/bin/time -v ./metagraph build -v \
            --parallel 30 -k 24 --mem-cap-gb 350 --suffix $F -o <GRAPH_DIR>/graph \
        2>&1 | tee <LOG_DIR>/log_$F"; \
done
```
2) Concatenate chunks
```bash
bsub -J StackChunks -W 12:00 -n 30 -R "rusage[mem=15000]" "/usr/bin/time -v \
    ~/metagraph concatenate -v -l 3 \
                                  -i <GRAPH_DIR>/graph \
                                  -o <GRAPH_DIR>/graph \
    2>&1 | tee <LOG_DIR>/log_stack.txt"
```

### Annotate graph
```bash
./metagraph annotate -v --anno-type row --fasta-anno \
                           -i primates.dbg \
                           -o primates \
                           ~/fasta_zurich/refs_chimpanzee_primates.fa
```

### Convert annotation to Multi-BRWT
1) Cluster columns
```bash
./metagraph transform_anno -v --linkage --greedy \
                           -o linkage.txt \
                           --subsample R \
                           -p NCORES \
                           primates.column.annodbg
```
Requires `N*R/8 + 6*N^2` bytes of RAM, where `N` is the number of columns and `R` is the number of rows subsampled.

2) Construct Multi-BRWT
```bash
./metagraph transform_anno -v -p NCORES --anno-type brwt \
                           -i linkage.txt \
                           -o primates \
                           --parallel-nodes V \
                           -p NCORES \
                           primates.column.annodbg
```
Requires `M*V/8 + Size(BRWT)` bytes of RAM, where `M` is the number of rows in the annotation and `V` is the number of nodes merged concurrently.

### Query graph
```bash
./metagraph query -v -i <GRAPH_DIR>/graph.dbg \
                        -a <GRAPH_DIR>/annotation.column.annodbg \
                        --discovery-fraction 0.8 --labels-delimiter ", " \
                        query_seq.fa
```

### Align to graph
```bash
./metagraph align -v -i <GRAPH_DIR>/graph.dbg query_seq.fa
```

### Assemble sequences
```bash
./metagraph assemble -v <GRAPH_DIR>/graph.dbg \
                        -o assembled.fa \
                        --unitigs
```

### Assemble differential sequences
```bash
./metagraph assemble -v <GRAPH_DIR>/graph.dbg \
                        --unitigs \
                        -a <GRAPH_DIR>/annotation.column.annodbg \
                        --label-mask-in LABEL_1 \
                        --label-mask-in LABEL_2 \
                        --label-mask-out LABEL_3 \
                        -o diff_assembled.fa
```

### Get stats
Stats for graph
```bash
./metagraph stats graph.dbg
```
Stats for annotation
```bash
./metagraph stats -a annotation.column.annodbg
```
Stats for both
```bash
./metagraph stats -a annotation.column.annodbg graph.dbg
```

## License
Metagraph is distributed under the GPLv3 License (see LICENSE).
Please find further information in the AUTHORS and COPYRIGHTS files.
