# Metagenome Graph Project

## Install

### Prerequisites
- cmake 3.6.1
- C++14
- HTSlib
- folly (optional)

All can be installed with `brew` or `linuxbrew`.

### Compile
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions  
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

3. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
4. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))`
5. run unit tests `./unit_tests`

### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (`Debug` by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)
- `-DPYTHON_INTERFACE=ON` -- compile python interface (requires shared libraries, OFF by default)

### Compile in the amino acid mode
- `-DCMAKE_DBG_ALPHABET=Protein` -- use the alphabet of amino acids (`DNA` by default)


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
./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 10 -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz 2>&1 | tee <LOG_DIR>/log.txt
```

* #### Chunking build
1) Build chunks
```bash
for F in {\$,A,C,G,T,N}{\$,A,C,G,T,N}; do ./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 100 -o <GRAPH_DIR>/graph --suffix $F <DATA_DIR>/*.fasta.gz 2>&1 | tee <LOG_DIR>/log_$F.txt; done
```
2) Concatenate chunks
```bash
./metagengraph concatenate -l 2 -i <GRAPH_DIR>/graph -o <GRAPH_DIR>/graph
```

#### Build from filtered reads
1) Filter reads
```bash
./metagengraph filter -v --parallel 30 -k 20 --noise-freq 3 <DATA_DIR>/*.fasta.gz 2>&1 | tee <LOG_DIR>/log_filter.txt
```
2) Build graph
```bash
./metagengraph build -v --parallel 30 -k 20 --mem-cap-gb 100 --noise-freq 3 -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz 2>&1 | tee <LOG_DIR>/log_build.txt
```

#### Distributed build
1) Build chunks
```bash
for F in {\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}{\\\$,A,C,G,T,N}; do bsub -J assemble$F -W 8:00 -n 30 -R "rusage[mem=15000]" "ls -1a <DATA_DIR>/*.fasta.gz | /usr/bin/time -v ./metagengraph build -v --parallel 30 -k 24 --mem-cap-gb 350 --suffix $F -o <GRAPH_DIR>/graph 2>&1 | tee <LOG_DIR>/log_$F"; done
```
2) Concatenate chunks
```bash
bsub -J StackChunks -W 12:00 -n 30 -R "rusage[mem=15000]" "/usr/bin/time -v ~/metagengraph concatenate -v -l 3 -i <GRAPH_DIR>/graph -o <GRAPH_DIR>/graph 2>&1 | tee <LOG_DIR>/log_stack.txt"
```

### Annotate graph
```bash
./metagengraph annotate -v --row-annotator --fasta-anno -i primates.dbg ~/fasta_zurich/refs_chimpanzee_primates.fa
```

### Query graph
```bash
./metagengraph classify -v -i <GRAPH_DIR>/graph.dbg -a <GRAPH_DIR>/annotation.color.annodbg --discovery-fraction 0.8 --labels-delimiter ", " query_seq.fa
```

### Get stats
Stats for graph
```bash
./metagengraph stats -a graph.dbg
```
Stats for annotation
```bash
./metagengraph stats -a annotation.color.annodbg
```
Stats for both
```bash
./metagengraph stats -a annotation.color.annodbg graph.dbg
```
