# Metagenome Graph Project

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ratschlab/metagraph)](https://github.com/ratschlab/metagraph/releases)
[![bioconda downloads](https://img.shields.io/conda/dn/bioconda/metagraph?color=blue)](https://bioconda.github.io/recipes/metagraph/README.html)
[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](#conda)
[![install with docker](https://img.shields.io/badge/install%20with-docker-brightgreen)](#docker)
[![install from source](https://img.shields.io/badge/install%20from-source-lightgrey)](#install-from-sources)
[![documentation](https://img.shields.io/badge/-online%20docs-grey)](https://metagraph.ethz.ch/static/docs/index.html)

MetaGraph is a tool for scalable construction of annotated genome graphs and sequence-to-graph alignment.

The default index representations in MetaGraph are extremely scalable and support building graphs with trillions of nodes and millions of annotation labels.
At the same time, the provided workflows and their careful implementation, combined with low-level optimizations of the core data structures, enable exceptional query and alignment performance.

#### Main features:
* Large-scale indexing of sequences
* [Python API](https://metagraph.ethz.ch/static/docs/api.html) for querying in the server mode
* Encoding [**k-mer counts**](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-counts) (e.g., expression values) and [**k-mer coordinates**](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-coordinates) (positions in source sequences)
* **Sequence alignment** against very large annotated graphs
* Scalable cleaning of very large de Bruijn graphs (to remove sequencing errors)
* Support for custom alphabets (e.g., {A,C,G,T,N} or amino acids)
* Algorithms for [differential assembly](https://metagraph.ethz.ch/static/docs/sequence_assembly.html#differential-assembly)

#### Design choices in MetaGraph:
* Use of succinct data structures and efficient representation schemes for extremely high scalability
* Algorithmic choices that work efficiently with succinct data structures (e.g., always prefer batched operations)
* Modular support of different graph and annotation representations
* Use of generic and extensible interfaces to support adding custom index representations / algorithms with little code overhead.

## Documentation
Online documentation is available at https://metagraph.ethz.ch/static/docs/index.html. Offline sources are [here](metagraph/docs/source).

## Install

### Conda

Install the [latest release](https://github.com/ratschlab/metagraph/releases/latest) on Linux or Mac OS X with Anaconda:

```
conda install -c bioconda -c conda-forge metagraph
```

### Docker

If docker is available on the system, immediately get started with

```
docker pull ghcr.io/ratschlab/metagraph:master
docker run -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master \
    build -v -k 10 -o /mnt/transcripts_1000 /mnt/transcripts_1000.fa
```
and replace `${HOME}` with a directory on the host system to map it under `/mnt` in the container.

To run the binary compiled for the `Protein` alphabet, just add `--entrypoint metagraph_Protein`:
```
docker run -v ${HOME}:/mnt --entrypoint metagraph_Protein ghcr.io/ratschlab/metagraph:master \
    build -v -k 10 -o /mnt/graph /mnt/protein.fa
```

As you see, running MetaGraph from docker containers is very easy. Also, the following command (or similar) may be handy to see what directory is mounted in the container or other sort of debugging of the command:
```
docker run -v ${HOME}:/mnt --entrypoint ls ghcr.io/ratschlab/metagraph:master /mnt
```

All different versions of the container image are listed [here](https://github.com/ratschlab/metagraph/pkgs/container/metagraph).

### Install From Sources

To compile from source (e.g., for builds with custom alphabet or other configurations), see [documentation online](https://metagraph.ethz.ch/static/docs/installation.html#install-from-source).


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

#### Build from k-mers filtered with KMC
```bash
K=20
./KMC/kmc -ci5 -t4 -k$K -m5 -fm <FILE>.fasta.gz <FILE>.cutoff_5 ./KMC
./metagraph build -v -p 4 -k $K --mem-cap-gb 10 -o graph <FILE>.cutoff_5.kmc_pre
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
                           --linkage-file linkage.txt \
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
                        --diff-assembly-rules diff_assembly_rules.json \
                        -o diff_assembled.fa
```

See [`metagraph/tests/data/example.diff.json`](metagraph/tests/data/example.diff.json) and [`metagraph/tests/data/example_simple.diff.json`](metagraph/tests/data/example_simple.diff.json) for sample files.

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

## Developer Notes

### Makefile

The `Makefile` in the top level source directory can be used to build and test `metagraph` more conveniently. The following
arguments are supported:
* `env`: environment in which to compile/run (`""`: on the host, `docker`: in a docker container)
* `alphabet`: compile metagraph for a certain alphabet (e.g. `DNA` or `Protein`, default `DNA`)
* `additional_cmake_args`: additional arguments to pass to cmake.

Examples:

```
# compiles metagraph in a docker container for the `DNA` alphabet
make build-metagraph env=docker alphabet=DNA
```


## License
Metagraph is distributed under the GPLv3 License (see LICENSE).
Please find further information in the AUTHORS and COPYRIGHTS files.
