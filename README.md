<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="metagraph/docs/source/images/metagraph_logo_dark.svg">
    <img src="metagraph/docs/source/images/metagraph_logo.svg" alt="MetaGraph" width="260">
  </picture>
</p>

# MetaGraph: Metagenome Graph Project

[![Platform: Linux | macOS](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-brightgreen)](#quick-start)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ratschlab/metagraph)](https://github.com/ratschlab/metagraph/releases)
[![Bioconda version](https://img.shields.io/conda/vn/bioconda/metagraph)](https://bioconda.github.io/recipes/metagraph/README.html)
[![bioconda downloads](https://img.shields.io/conda/dn/bioconda/metagraph?color=blue)](https://bioconda.github.io/recipes/metagraph/README.html)
[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](https://metagraph.ethz.ch/static/docs/installation.html#install-with-conda)
[![install with docker](https://img.shields.io/badge/install%20with-docker-brightgreen)](#docker)
[![install from source](https://img.shields.io/badge/install%20from-source-lightgrey)](#install-from-source)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41586--025--09603--w-blue)](https://doi.org/10.1038/s41586-025-09603-w)
[![documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://metagraph.ethz.ch/static/docs/index.html)

**Scalable indexing and querying of annotated genome graphs, from a handful of genomes to petabase-scale sequence repositories.**

Index your reads or assemblies once. Then search them for matching samples, alignments and source positions.

The MetaSUB index on [MetaGraph Online](https://metagraph.ethz.ch/search) covers 7.2 trillion bases from about 4,200 samples. In a direct-server test, it returned matches for a 10 kbp DNA query in under a second, using four query threads on a warm index.

A MetaGraph index has 2 parts. The **de Bruijn graph** stores every *k*-mer of the input sequences. The **annotation matrix** links each *k*-mer to its labels, such as the sample, the sequence header, a count or a position.

```mermaid
flowchart LR
    F["FASTA / FASTQ<br/>reads, assemblies, transcripts"]
    F -->|build| G[("de Bruijn graph<br/>graph.dbg")]
    F -->|annotate| M[("annotation matrix<br/>*.annodbg")]
    G --> Q(("query / align"))
    M --> Q
    Qseq["query sequence"] --> Q
    Q --> R["matches<br/>labels, counts, positions"]
```

### Features

- **Search public archives** at [metagraph.ethz.ch](https://metagraph.ethz.ch). See [MetaGraph Online](#metagraph-online).
- **Index your own data.** Build a *k*-mer index over reads, assemblies or transcripts. Then query it for matching labels.
- **Optional per-*k*-mer payloads.** Attach [counts](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-counts) (abundance, e.g. expression or coverage) or [coordinates](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-coordinates) (positions that losslessly encode source sequences and return per-target hit positions).
- **Sequence alignment.** Align sequences to the full annotated graph, with sub-*k* seeding for arbitrarily short queries.
- **Scalable graph cleaning.** Strip sequencing errors out of very large de Bruijn graphs.
- **[Differential assembly](https://metagraph.ethz.ch/static/docs/sequence_assembly.html#differential-assembly).** Extract sequences present in one group of samples and absent from another.
- **[Python API & HTTP server](https://metagraph.ethz.ch/static/docs/api.html).** Drive MetaGraph from Python or query a running instance over HTTP.

<details>
<summary>Under the hood</summary>

- **Succinct data structures**: the default `succinct` (BOSS) graph representation uses only 2–4 bits per *k*-mer.
- **Modular annotation formats**: `ColumnCompressed`, `RowDiff<Multi-BRWT>`, `RowSparse`, `Rainbowfish`, plus count- and coordinate-aware variants. Pick the compression/speed tradeoff that fits your scale.
- **Custom alphabets**: `{A,C,G,T}`, `{A,C,G,T,N}`, amino acids, case-sensitive DNA, or compile-time custom alphabets.
- **Extensible by design**: generic interfaces let developers add new graph/annotation representations or algorithms with little code.
- **Memory-mapped loading**: pass `--mmap` to any command to start fast and use less RAM at query time. NVMe works best. An SSD works but is slower.
- **Scales to trillions of *k*-mers and millions of labels**: the hosted indexes at metagraph.ethz.ch cover all public SRA data.

</details>

> **Full documentation:** <https://metagraph.ethz.ch/static/docs/index.html> (offline copy in [`metagraph/docs/source`](metagraph/docs/source)).

## MetaGraph Online

Try MetaGraph without installing anything at <https://metagraph.ethz.ch/search>. It searches indexes of public DNA, RNA and protein data, such as RefSeq, UHGG, Tara Oceans and UniParc. See the [list of databases](https://metagraph.ethz.ch/indexes). Paste a sequence and pick the indexes to search.

Prefer local analysis? The prebuilt indexes are also published on [AWS Open Data](https://registry.opendata.aws/metagraph/) for download (no AWS account needed) and offline querying. See [Preconstructed indexes](https://metagraph.ethz.ch/static/docs/resources.html#preconstructed-indexes).

## Quick start

### 1. Start the example environment

The Docker image includes the MetaGraph CLI, the workflow wrapper and the example FASTA files.
From a directory where you want to keep the index, start a shell in the image:

```bash
mkdir -p out
docker run --rm -it -v "$PWD/out:/opt/metagraph/out" -w /opt/metagraph ghcr.io/ratschlab/metagraph:master bash
```

Run the following build and query commands in that shell. The index stays in the host's `out/` directory after you exit.

### 2. Build an index

```bash
metagraph-workflows build <(ls metagraph/tests/data/*.fa) -o out/ --primary -p 2 --mem-gb 1
```

`--primary` indexes 1 strand of each *k*-mer pair. The index is about half the size. Use it when the strand of the reads is unknown, as in most short-read sequencing.

This runs the full pipeline. It writes the graph `graph.dbg`, a smaller and slower copy `graph_small.dbg`, and the annotation `graph.relax.row_diff_brwt.annodbg`. The [Quick start guide](https://metagraph.ethz.ch/static/docs/quick_start.html) explains each step.

<details>
<summary>Real-workload example with file list and hardware budget</summary>

```bash
metagraph-workflows build samples.txt -o out/ --primary \
  -p 34 --mem-gb 70 --disk-swap-dir /scratch/swap
```

The positional argument accepts a file list (one path per line), a directory, or process substitution. Per-stage memory caps, thread packing, and BRWT parameters are derived from the budget. See `metagraph-workflows build --help` for all options.

</details>

### 3. Query

Query the index. Any FASTA file works as input. The bundled test data is enough for a first try.

```bash
metagraph query --query-mode matches -p 2 \
    -i out/graph.dbg -a out/graph.relax.row_diff_brwt.annodbg \
    metagraph/tests/data/transcripts_100.fa
```

Other ways to use the index:

- [`metagraph align`](https://metagraph.ethz.ch/static/docs/sequence_search.html#sequence-to-graph-alignment): align sequences to the graph and report the alignment itself. With a coordinate-aware annotation it acts as a read mapper, returning source positions.
- `metagraph query --align`: align each query to the graph first, then report the *labels* of the best-scoring alignment instead of requiring exact *k*-mer matches. Useful for divergent or noisy queries.
- [`metagraph server_query`](https://metagraph.ethz.ch/static/docs/api.html): start an HTTP server that answers queries over a REST API (also drivable from the Python client).

The [Minimal example](#minimal-example) below walks through each step on a smaller dataset.

## Minimal example

<details>
<summary>Step-by-step walkthrough with the <code>metagraph</code> CLI (build → annotate → query → stats)</summary>

A hands-on demo using `metagraph` directly (no workflow wrapper). A *label* is whatever tag you want each *k*-mer associated with: a filename, a FASTA header, or a custom string. `--anno-header` below produces one label per FASTA record.

```bash
cd metagraph/tests/data

# 1. Build a de Bruijn graph (k-mer index) with k=31
metagraph build -v -p 4 -k 31 -o samples metasub_fake_data_simple.fa

# 2. Construct an annotation matrix linking each k-mer to its source label.
metagraph annotate -v -p 4 -i samples.dbg --anno-header \
                   -o samples metasub_fake_data_simple.fa

# 3. Query: for each query sequence, report all labels whose k-mers cover ≥80% of it.
metagraph query --query-mode matches \
                -i samples.dbg -a samples.column.annodbg \
                --min-kmers-fraction-label 0.8 \
                metasub_fake_data_simple.fa

# 4. Print graph + annotation stats.
metagraph stats -a samples.column.annodbg samples.dbg
```

Outputs `samples.dbg` (the de Bruijn graph) and `samples.column.annodbg` (3 labels, one per FASTA record). Sample query output (tab-separated: query index, query header, then one `<label>:k-mer-count` entry per matching label):

```
0   kl_sample   <kl_sample>:330
1   zh_sample   <kl_sample>:243   <zh_sample>:243
2   tk_sample   <kl_sample>:207   <tk_sample>:207
```

Each query matches at least its own label. `zh_sample` has 243 *k*-mers, and all of them are in `kl_sample`. The same holds for `tk_sample` (207 of 207). Both clear the 80% threshold. `kl_sample` matches only itself, because neither shorter sample covers its 330 *k*-mers.

</details>

## More recipes

<details>
<summary>Direct CLI: align / assemble / differential assembly / stats</summary>

```bash
# Build a de Bruijn graph from a fasta (single sample)
metagraph build -v -p 8 -k 31 --mem-cap-gb 10 -o graph data.fa.gz

# Build using disk swap (limits RAM at the cost of disk I/O)
metagraph build -v -p 8 -k 31 --mem-cap-gb 10 --disk-swap /scratch/swap -o graph data.fa.gz

# Annotate a graph with file-based labels (one label per input file)
metagraph annotate -v -p 8 --anno-filename -i graph.dbg -o annotation data.fa.gz

# Align sequences to the graph (plain sequence-to-graph alignment, no labels).
metagraph align -v -i graph.dbg query.fa

# Labeled alignment: with -a, the walk continues only while its k-mers share
# at least 1 label. It reports those shared labels.
metagraph align -v -i graph.dbg -a annotation.row_diff_brwt.annodbg reads.fq

# Same, but with a coordinate-aware annotator: the walk is additionally
# coordinate-consistent (positions step by ±1 per node), so hits map to
# source positions. Functions as a read mapper over indexed genomes.
metagraph align -v -i graph.dbg -a annotation.row_diff_brwt_coord.annodbg reads.fq

# query --align: aligns each query to the graph WITHOUT label constraints,
# then fetches labels for the highest-scoring walk. Use when exact k-mer
# matching is too strict (divergent or noisy queries).
metagraph query --align -i graph.dbg -a annotation.row_diff_brwt.annodbg query.fa

# Assemble unitigs from a graph (writes assembled.fasta.gz)
metagraph assemble -v graph.dbg -o assembled --unitigs

# Differential assembly: JSON rules define which label groups must be present vs. absent.
# Sample rule files: metagraph/tests/data/example.diff.json, example_simple.diff.json.
metagraph assemble -v graph.dbg --unitigs \
    -a annotation.column.annodbg \
    --diff-assembly-rules diff_assembly_rules.json \
    -o diff_assembled

# Stats: graph only, annotation only, or both
metagraph stats graph.dbg
metagraph stats -a annotation.column.annodbg
metagraph stats -a annotation.column.annodbg graph.dbg
```

See the [online docs](https://metagraph.ethz.ch/static/docs/index.html) for the full subcommand reference, alphabet builds, and the row-diff / BRWT clustering pipeline.

</details>

<details>
<summary>Build from KMC k-mer counters (e.g., to filter low-abundance k-mers)</summary>

[KMC](https://github.com/refresh-bio/KMC) is a fast *k*-mer counter that can pre-process inputs and filter by abundance.

Configure with `cmake .. -DBUILD_KMC=ON` to build the standalone KMC executable
used below. Reading existing KMC files is always supported.

Build a graph from *k*-mers occurring at least 5 times:

```bash
K=31
# Count k-mers with KMC (drops k-mers with abundance < 5)
./KMC/kmc -ci5 -t4 -k$K -m5 -fq input.fastq.gz input.cutoff_5 ./KMC

# Build the graph directly from the KMC counter
metagraph build -v -p 4 -k $K -o graph input.cutoff_5.kmc_pre
```

Add `--count-kmers` to keep the KMC abundance counts as a weight vector alongside the graph (`graph.dbg.weights`). The weighted graph is the input for [`metagraph clean`](https://metagraph.ethz.ch/static/docs/quick_start.html#graph-cleaning) and for indexing *k*-mer counts.

</details>

<details>
<summary>Common query flags</summary>

```bash
# Filter labels with low k-mer coverage (default 0.7)
metagraph query --min-kmers-fraction-label 0.8 ...

# Delimiter joining labels in `--query-mode labels` output (default ":")
metagraph query --query-mode labels --labels-delimiter ", " ...

# Per-k-mer presence/absence bitmask per matching label
metagraph query --query-mode signature ...

# Indexed k-mer counts (requires count-aware annotation, see online docs)
metagraph query --query-mode counts ...

# k-mer coordinates / per-target positions (requires coordinate-aware annotation)
metagraph query --query-mode coords ...

# Server mode (Python API / HTTP queries)
metagraph server_query -i graph.dbg -a annotation.row_diff_brwt.annodbg --port 5555 --parallel 8
```

</details>

## More install options

MetaGraph runs on Linux and macOS. Bioconda ships the `DNA` and `Protein` alphabets. `metagraph` is a link to `metagraph_DNA`. The Docker image adds `DNA5`. For other alphabets, build from source.

### Docker

```bash
docker pull ghcr.io/ratschlab/metagraph:master
docker run -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master \
    metagraph build -v -k 10 -o /mnt/transcripts_1000 /mnt/transcripts_1000.fa
```

Replace `${HOME}` with the host directory you want to expose under `/mnt`. The default image uses the `DNA` alphabet `{A,C,G,T}`. The same image runs `metagraph_DNA5` and `metagraph_Protein`. All published images are listed [here](https://github.com/ratschlab/metagraph/pkgs/container/metagraph).

<details>
<summary>Advanced Docker usage</summary>

Inspect what's mounted in the container:

```bash
docker run -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master ls /mnt
```

Run the `DNA5` / `Protein` variants:

```bash
docker run -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master \
    metagraph_Protein build -v -k 10 -o /mnt/graph /mnt/protein.fa
```

Drop into an interactive shell for multi-step workflows:

```bash
docker run -it --entrypoint /bin/bash -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master
# inside container:
metagraph --version
metagraph build -v -k 31 -o /mnt/graph /mnt/data.fa.gz
```

</details>

### Install from source

See the [installation guide](https://metagraph.ethz.ch/static/docs/installation.html#install-from-source) for custom alphabet / debug builds.

## Contributing

<details>
<summary>Developer notes: docker build, releases</summary>

**Build a docker container.** Run `docker build .`

**Releases.** 1. Bump the version in `package.json`. 2. Tag the commit with that version. 3. Create a GitHub release.

</details>

## Citation

If MetaGraph or its index resources are useful in your work, please cite:

> Karasikov M, Mustafa H, Danciu D, Kulkov O, Zimmermann M, Barber C, Rätsch G, Kahles A. Efficient and accurate search in petabase-scale sequence repositories. *Nature*. 2025;647: 1036–1044. <https://www.nature.com/articles/s41586-025-09603-w>

<details>
<summary>BibTeX</summary>

```bibtex
@article{karasikov2025metagraph,
  title={Efficient and accurate search in petabase-scale sequence repositories},
  author={Karasikov, Mikhail and Mustafa, Harun and Danciu, Daniel and Kulkov, Oleksandr and Zimmermann, Marc and Barber, Christopher and R{\"a}tsch, Gunnar and Kahles, Andr{\'e}},
  journal={Nature},
  volume={647},
  number={8091},
  pages={1036--1044},
  year={2025},
  publisher={Nature Publishing Group},
  doi={10.1038/s41586-025-09603-w}
}
```

</details>

## License

MetaGraph is distributed under the GPLv3 License (see [LICENSE](LICENSE)). See also [AUTHORS](AUTHORS) and [COPYRIGHT](COPYRIGHT).
