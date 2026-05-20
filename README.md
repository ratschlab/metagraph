<p align="center">
  <img src="metagraph/docs/source/images/metagraph_logo.png" alt="MetaGraph" width="420">
</p>

[![Linux](https://img.shields.io/badge/Linux-supported-brightgreen?logo=linux&logoColor=white)](#-quick-start)
[![macOS](https://img.shields.io/badge/macOS-supported-brightgreen?logo=apple&logoColor=white)](#-quick-start)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ratschlab/metagraph)](https://github.com/ratschlab/metagraph/releases)
[![Bioconda version](https://img.shields.io/conda/vn/bioconda/metagraph)](https://bioconda.github.io/recipes/metagraph/README.html)
[![bioconda downloads](https://img.shields.io/conda/dn/bioconda/metagraph?color=blue)](https://bioconda.github.io/recipes/metagraph/README.html)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41586--025--09603--w-blue)](https://doi.org/10.1038/s41586-025-09603-w)
[![documentation](https://img.shields.io/badge/📖-online%20docs-blue.svg)](https://metagraph.ethz.ch/static/docs/index.html)

**Scalable indexing and querying of annotated genome graphs — from a handful of genomes to petabase-scale sequence repositories.**

MetaGraph builds compact, searchable indexes over sequencing data: assemblies, raw reads, transcripts, whole bacterial collections. Once indexed, you can query for sequences in milliseconds, align reads against trillions of k-mers, and recover where every match came from.

### Features

- ⚡ **Scale.** Indexes with trillions of k-mers and millions of annotation labels; petabase-scale collections of public sequencing data have been indexed end-to-end.
- 🔢 **[k-mer counts](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-counts)** and 📏 **[k-mer coordinates](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-coordinates).** Optional payloads alongside k-mer presence — for expression levels, alignment-position recovery, or losslessly encoding source genomes.
- 🧬 **Sequence alignment** against the full annotated graph, with sub-k seeding for arbitrarily short queries.
- 🧹 **Scalable graph cleaning** to strip sequencing errors out of very large de Bruijn graphs.
- 🔤 **Custom alphabets** — `{A,C,G,T}`, `{A,C,G,T,N}`, protein, or your own.
- 🔀 **[Differential assembly](https://metagraph.ethz.ch/static/docs/sequence_assembly.html#differential-assembly)** — extract sequences present in one group of samples and absent from another.

> 📖 **Full documentation:** <https://metagraph.ethz.ch/static/docs/index.html>

## 🌐 MetaGraph Online

Try MetaGraph without installing anything: <https://metagraph.ethz.ch> hosts a public search engine over 50+ petabases of DNA, RNA, and protein sequences, indexing SRA, ENA, DRA, RefSeq, UniProt, and more.

## 🚀 Quick start

### 1. Install

```bash
conda create -n metagraph python
conda activate metagraph
conda install -c bioconda -c conda-forge metagraph
pip install --force-reinstall "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/workflows"
```

### 2. Build an index

Clone the repo for the bundled test data, then build:

```bash
git clone https://github.com/ratschlab/metagraph.git && cd metagraph
metagraph-workflows build <(ls metagraph/tests/data/*.fa) -o out/ --primary
```

The positional argument accepts a file list (one path per line), a directory, or — as above — process substitution.

For real workloads, supply your own sample list and a hardware budget:

```bash
metagraph-workflows build samples.txt -o out/ --primary \
  -p 34 --mem-gb 70 --disk-swap-dir /scratch/swap
```

Per-stage memory caps, thread packing, and BRWT parameters are derived from the budget. See `metagraph-workflows build --help` for all options.

Outputs: `graph.dbg`, `graph_small.dbg`, and one `.annodbg` per requested format.

Under the hood, the workflow runs four `metagraph` stages — see [the pipeline docs](https://metagraph.ethz.ch/static/docs/quick_start.html) for each:

1. `metagraph build` — graph construction
2. `metagraph annotate` — per-sample annotation
3. `metagraph transform_anno --anno-type row_diff` — row-diff transform (stages 0–2)
4. `metagraph transform_anno --anno-type *_brwt` — BRWT clustering

### 3. Query

The default `--anno-type` produces `<graph>.relax.row_diff_brwt.annodbg`. Query it against the graph (any fasta works as input — the bundled test data is fine for a smoke test):

```bash
metagraph query --query-mode matches -p 8 \
    -i out/graph.dbg -a out/graph.relax.row_diff_brwt.annodbg \
    metagraph/tests/data/transcripts_100.fa
```

For the [Python API](https://metagraph.ethz.ch/static/docs/api.html), see the online docs.

<details>
<summary>🔧 Direct <code>metagraph</code> CLI (align / assemble / stats / ...)</summary>

The workflow wrapper covers the common build path. If you'd rather invoke `metagraph` directly:

```bash
# Build a graph from a fasta (single sample, simplest case)
metagraph build -v -p 8 -k 31 --mem-cap-gb 10 -o graph data.fa.gz

# Build using disk swap (limits RAM at the cost of disk I/O)
metagraph build -v -p 8 -k 31 --mem-cap-gb 10 --disk-swap /scratch/swap -o graph data.fa.gz

# Annotate a graph with file-based labels
metagraph annotate -v -p 8 --anno-filename -i graph.dbg -o annotation data.fa.gz

# Align reads against a graph
metagraph align -v -i graph.dbg query.fa

# Assemble unitigs from a graph
metagraph assemble -v graph.dbg -o assembled.fa --unitigs

# Differential assembly (extract sequences present in some labels, absent in others)
metagraph assemble -v graph.dbg --unitigs \
    -a annotation.column.annodbg \
    --diff-assembly-rules diff_assembly_rules.json \
    -o diff_assembled.fa

# Stats for graph + annotation
metagraph stats -a annotation.column.annodbg graph.dbg
```

See the [online docs](https://metagraph.ethz.ch/static/docs/index.html) for the full subcommand reference, alphabet builds, and the row-diff / BRWT clustering pipeline.

</details>

## 📦 Install

The recommended conda + pip install is covered in [Quick start](#-quick-start). Alternative install methods:

### 🐳 Docker

```bash
docker pull ghcr.io/ratschlab/metagraph:master
docker run -v ${HOME}:/mnt ghcr.io/ratschlab/metagraph:master \
    metagraph build -v -k 10 -o /mnt/transcripts_1000 /mnt/transcripts_1000.fa
```

Replace `${HOME}` with the host directory you want to expose under `/mnt`. The default image targets the `DNA` alphabet `{A,C,G,T}`; the `DNA5` and `Protein` variants are invoked as `metagraph_DNA5` / `metagraph_Protein`. All published images are listed [here](https://github.com/ratschlab/metagraph/pkgs/container/metagraph).

### 🛠️ Install from sources

See the [installation guide](https://metagraph.ethz.ch/static/docs/installation.html#install-from-source) for custom alphabet / debug builds.

## 📝 Citation

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

## ⚖️ License

MetaGraph is distributed under the GPLv3 License (see [LICENSE](LICENSE)). See also [AUTHORS](AUTHORS) and [COPYRIGHTS](COPYRIGHTS).
