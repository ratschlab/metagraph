# Metagenome Graph Project

## Install

### Prerequisites
- cmake 3.6.1
- GNU GCC (C++14)
- HTSlib
- folly (optional)

All can be installed with `brew` or `linuxbrew`.

### Compile
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. make sure all submodules are downloaded: `git submodule update --init --recursive`
3. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions  
or simply run the following script
```bash
git submodule update --init --recursive
cd metagraph/external-libraries/sdsl-lite
./install.sh $(pwd)

cd ../libmaus2
libtoolize
aclocal
autoreconf -i -f
./configure --prefix=$(pwd)
make -j $(($(getconf _NPROCESSORS_ONLN) - 1))
make install
cd ../../../
```
use `glibtoolize` instead of `libtoolize` on MacOS

4. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
5. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))`
6. run unit tests `./unit_tests`

#### Typical issues
* Linking against dynamic libraries in Anaconda when compiling libmaus2 (make sure that packages like Anaconda are not listed in the exported environment variables).

### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (`Debug` by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)
- `-DPYTHON_INTERFACE=ON` -- compile python interface (requires shared libraries, OFF by default)
- `-DBUILD_KMC=OFF` -- do not compile the KMC executable (ON by default)
- `-DWITH_AVX=OFF` -- compile without support for the avx instructions (ON by default)
- `-DCMAKE_DBG_ALPHABET=[Protein|DNA|DNA_CASE_SENSITIVE]` -- alphabet to use (`DNA` by default)


## Print usage
```bash
./metagengraph
```
```bash
./metagengraph build
```

## Usage examples

For real examples, see [scripts](./scripts).

### Build graph

* #### Simple build
```bash
./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 10 \
                        -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz \
2>&1 | tee <LOG_DIR>/log.txt
```

* #### Chunking build
1) Build chunks
```bash
for F in {\$,A,C,G,T,N}{\$,A,C,G,T,N}; do \
    ./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 100 \
                            -o <GRAPH_DIR>/graph --suffix $F \
                            <DATA_DIR>/*.fasta.gz \
    2>&1 | tee <LOG_DIR>/log_$F.txt; \
done
```
2) Concatenate chunks
```bash
./metagengraph concatenate -l 2 -i <GRAPH_DIR>/graph -o <GRAPH_DIR>/graph
```

#### Build from filtered reads
1) Filter reads
  * using filtering in blocks
```bash
./metagengraph filter -v --parallel 30 -k 20 --filter-abund 3 <DATA_DIR>/*.fasta.gz
```
  * using KMC
```bash
./KMC/kmc -k21 -m5 -fq -t30 <FILE>.fasta.gz <FILE>.fasta.gz.kmc ./KMC
./metagengraph filter -v --parallel 30 -k 20 --filter-abund 3 --kmc <FILE>.fasta.gz
```
2) Build graph
```bash
./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 100 --filter-abund 3 \
                        -o <GRAPH_DIR>/graph \
                        <DATA_DIR>/*.fasta.gz
```

#### Distributed build
1) Build chunks
```bash
for F in {\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}; do \
    bsub -J assemble$F -W 8:00 -n 30 -R "rusage[mem=15000]" \
        "ls -1a <DATA_DIR>/*.fasta.gz | /usr/bin/time -v ./metagengraph build -v \
            --parallel 30 -k 24 --mem-cap-gb 350 --suffix $F -o <GRAPH_DIR>/graph \
        2>&1 | tee <LOG_DIR>/log_$F"; \
done
```
2) Concatenate chunks
```bash
bsub -J StackChunks -W 12:00 -n 30 -R "rusage[mem=15000]" "/usr/bin/time -v \
    ~/metagengraph concatenate -v -l 3 \
                                  -i <GRAPH_DIR>/graph \
                                  -o <GRAPH_DIR>/graph \
    2>&1 | tee <LOG_DIR>/log_stack.txt"
```

### Annotate graph
```bash
./metagengraph annotate -v --row-annotator --fasta-anno \
                           -i primates.dbg \
                           ~/fasta_zurich/refs_chimpanzee_primates.fa
```

### Query graph
```bash
./metagengraph classify -v -i <GRAPH_DIR>/graph.dbg \
                           -a <GRAPH_DIR>/annotation.column.annodbg \
                           --discovery-fraction 0.8 --labels-delimiter ", " \
                           query_seq.fa
```

### Get stats
Stats for graph
```bash
./metagengraph stats graph.dbg
```
Stats for annotation
```bash
./metagengraph stats -a annotation.column.annodbg
```
Stats for both
```bash
./metagengraph stats -a annotation.column.annodbg graph.dbg
```
