# pymetagraph

Load a MetaGraph index in Python for label queries and seed-and-extend alignment.
This package combines the generic query API with optional labeled alignment; it
has no Readfish dependency. For remote servers, use the separate `metagraph`
client package in `metagraph/api/python`.

## Build

Use the repository's C++20 toolchain and dependencies, including its initialized
git submodules. The current SDSL dependency is header-only and does not require a
separate install step.

```sh
python -m pip install ./metagraph/api/python/pymetagraph
```

The package enables `BUILD_PYTHON_BINDINGS` and builds position-independent libraries.
Linux wheels bundle shared Boost libraries in `pymetagraph/.libs`; other platforms
link Boost statically. The system C++ and compression-library runtimes are still
required. Wheels built here target the build environment; use tools such as
`auditwheel` when preparing portable Linux distributions.
Jemalloc is disabled in this configuration to avoid static
TLS failures when importing the extension. Other build directories keep their
normal allocator configuration. Building the native library can take several minutes.

For direct CMake development, use a Python environment with `pybind11` installed:

```sh
python -m pip install pybind11
cmake -S metagraph -B build-python -DBUILD_PYTHON_BINDINGS=ON -DMETAGRAPH_BUILD_EXECUTABLE=OFF -DPython_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build build-python --target _pymetagraph_core -j8
cmake --install build-python --component _pymetagraph_core --prefix "$PWD/metagraph/api/python/pymetagraph"
export PYTHONPATH="$PWD/metagraph/api/python/pymetagraph${PYTHONPATH:+:$PYTHONPATH}"
python -c "from pymetagraph import AnnotatedGraph"
```

The native module must be installed alongside the Python files in the `pymetagraph`
package. The pip command handles that layout automatically.

## Query

```python
from pymetagraph import AnnotatedGraph

graph = AnnotatedGraph("graph.dbg", "graph.column.annodbg", threads=4)
graph.k
graph.labels
graph.query("ACGT...", num_top_labels=5, discovery_fraction=0.1)
# [(label, matching_kmer_count), ...]
graph.query(["ACGT...", "TTGA..."])  # one result list per sequence
```

`query` has the same label-count semantics as `metagraph query`. A sequence shorter
than k or with no qualifying hits gets an empty list. Annotation is required for
this method.

## Align

```python
graph.align("ACGT...", only_forwards=True, xdrop=27, rel_score_cutoff=0.95)
# [{"score": ..., "strand": "+", "query_start": ..., "query_end": ...,
#   "annotation": {"cigar": ..., "ref_sequence": ..., "labels": [...]}, ...}]

unlabeled = AnnotatedGraph("graph.dbg", threads=4)
unlabeled.align(["ACGT...", "TTGA..."])
```

Alignment results use the CLI JSON fields, without graph-path mappings, plus
`strand` and `query_start`/`query_end`. Query positions are zero-based, half-open
intervals in the original input. Reference `nt_coords` inside label entries retain
the CLI's one-based inclusive format; an available `.seqs` index resolves them to
sequence headers. Query positions must not be used as reference coordinates.

A string input returns a list of alignments; a batch returns one such list per
input. No alignment gives an empty list. Both methods release the GIL while the
native work runs and use the configured worker pool.

Alignment options include `seed_length` (0 selects min(19, k)),
`max_alternative_alignments`, `min_exact_match`, `xdrop`,
`rel_score_cutoff`, `only_forwards`, `connect_anchors`,
`extend_chains`, `chain_alignments`, and `two_sided`. Chaining and two-sided
extension default to false. Enable both for the more sensitive search;
`only_forwards` independently selects the query orientations. Scoring defaults
match the CLI.
With `connect_anchors=False`, a result containing an unresolved graph gap raises
`ValueError`: raw chains cannot use the alignment JSON schema. Keep anchor
connection enabled to obtain inexact alignments. Batch failures drain all workers
before raising, and the loaded graph remains usable.

Search limits use the CLI defaults: `max_num_seeds_per_locus=1000`,
`max_num_nodes_per_seq_char=5`, and `max_ram_per_alignment=200` (decimal MB).
Overfull seed loci are discarded. Each connection or extension window stops at
its graph-position or allocated-workspace limit and retains valid partial results;
a limit can therefore reduce alignment sensitivity. The RAM cap covers the
dynamic-programming tables, hash buckets, search queue and visited-position set,
including allocated capacity. A graph position is a node and its distance from
the seed; repeated scoring visits count once. It excludes loaded indices,
seed/annotation caches and result storage,
so it is not a process RSS limit. Zero node/RAM limits disable extension search.

`rel_score_cutoff` filters completed alignments relative to the best score for
each label (or the best unlabeled result). Seed-chain scores are not used for
this filter because extension can improve them.

## Tests and limits

```sh
METAGRAPH=/path/to/metagraph python -m pytest metagraph/api/python/pymetagraph/tests
```

The tests construct their own indices and compare both query and alignment output
with the CLI. One annotation file is supported. Indices must be compatible with
the linked MetaGraph version: missing paths are rejected by the Python wrapper,
but the existing native loaders may terminate the process for corrupt or
incompatible files. Protein alignment and performance have not been qualified
by these DNA tests.
