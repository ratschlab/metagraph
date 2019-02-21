# Metagenome Graph Project

## Install

### Prerequisites
- cmake 3.6.1
- GNU GCC with C++17 (gcc-8 or higher) or LLVM Clang (clang-7 or higher)
- HTSlib
- jsoncpp
- boost
- folly (optional)

All can be installed with [brew](https://brew.sh) or [linuxbrew](https://linuxbrew.sh)

#### For compiling with GNU GCC:
```
brew install gcc autoconf automake libtool cmake make htslib jsoncpp
brew install --build-from-source boost
(optional) brew install --build-from-source double-conversion gflags glog lz4 snappy zstd folly
brew install gcc@8
```
Then set the environment variables accordingly:
```
echo "\
# Use gcc-8 with cmake
export CC=\"\$(which gcc-8)\"
export CXX=\"\$(which g++-8)\"
" >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )
```

#### For compiling with LLVM Clang:
```
brew install llvm libomp autoconf automake libtool cmake make htslib jsoncpp boost folly
```
Then set the environment variables accordingly:
```
echo "\
# OpenMP
export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix libomp)/lib\"
export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix libomp)/include\"
export CXXFLAGS=\"\$CXXFLAGS -I$(brew --prefix libomp)/include\"
# Clang C++ flags
export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix llvm)/lib -Wl,-rpath,$(brew --prefix llvm)/lib\"
export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix llvm)/include\"
export CXXFLAGS=\"\$CXXFLAGS -stdlib=libc++\"
# Path to Clang
export PATH=\"$(brew --prefix llvm)/bin:\$PATH\"
# Use Clang with cmake
export CC=\"\$(which clang)\"
export CXX=\"\$(which clang++)\"
" >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )
```


### Compile
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. make sure all submodules are downloaded: `git submodule update --init --recursive`
3. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions  
or simply run the following script
```bash
git submodule update --init --recursive

pushd metagraph/external-libraries/sdsl-lite
./install.sh $(pwd)
popd

pushd metagraph/external-libraries/libmaus2
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(pwd) .
make -j $(($(getconf _NPROCESSORS_ONLN) - 1))
make install
popd
```

4. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
5. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))`
6. run unit tests `./unit_tests`

#### Typical issues
* Linking against dynamic libraries in Anaconda when compiling libmaus2
  * Solution: make sure that packages like Anaconda are not listed in the exported environment variables

* Trying to link against dynamic OpenMP libraries when compiling in static executables.
  * Solution: apply patch
```diff
--- a/metagraph/CMakeLists.txt
+++ b/metagraph/CMakeLists.txt
@@ -42,6 +42,10 @@ if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
 endif()

 find_package(OpenMP REQUIRED)
+if(OPENMP_FOUND)
+    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
+    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
+endif()
 if(NOT TARGET OpenMP::OpenMP_CXX)
     add_library(OpenMP_TARGET INTERFACE)
     add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
@@ -188,7 +192,6 @@ set(METALIBS ${METALIBS}
   -lssl -lcrypto -llzma
   -lbrwt
   -latomic
-  OpenMP::OpenMP_CXX
 )

 if(BUILD_STATIC)
 ```


### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (`Release` by default)
- `-DBUILD_STATIC=[ON|OFF]` -- link statically (`OFF` by default)
- `-DPYTHON_INTERFACE=[ON|OFF]` -- compile python interface (requires shared libraries, `OFF` by default)
- `-DBUILD_KMC=[ON|OFF]` -- compile the KMC executable (`ON` by default)
- `-DWITH_AVX=[ON|OFF]` -- compile with support for the avx instructions (`ON` by default)
- `-DCMAKE_DBG_ALPHABET=[Protein|DNA|DNA4|DNA_CASE_SENSITIVE]` -- alphabet to use (`DNA` by default)

## Typical workflow
1. Build de Bruijn graph from Fasta files, FastQ files, or [KMC k-mer counters](https://github.com/refresh-bio/KMC/):\
`./metagengraph build`
2. Annotate graph using the column compressed annotation:\
`./metagengraph annotate`
3. Transform the built annotation to a different annotation scheme:\
`./metagengraph transform_anno`
4. Merge annotations (optional):\
`./metagengraph merge_anno`
5. Query annotated graph\
`./metagengraph classify`

### Example
```
DATA="../tests/data/transcripts_1000.fa"

./metagengraph build -k 12 -o transcripts_1000 $DATA

./metagengraph annotate -i transcripts_1000 --anno-filename -o transcripts_1000 $DATA

./metagengraph classify -i transcripts_1000 -a transcripts_1000.column.annodbg $DATA

./metagengraph stats -a transcripts_1000 transcripts_1000
```

For real examples, see [scripts](./metagraph/scripts).

### Print usage
`./metagengraph`

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

#### Build with filtering k-mers using KMC
```bash
CUTOFF=5
K=20
./KMC/kmc -ci$CUTOFF -t30 -k$K -m5 -fq -b <FILE>.fasta.gz <FILE>.kmc_$CUTOFF ./KMC
./metagengraph build -v -p 30 -k $K --mem-cap-gb 10 --kmc -o graph <FILE>.kmc_$CUTOFF
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
