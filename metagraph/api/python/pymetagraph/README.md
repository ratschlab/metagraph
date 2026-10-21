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

The package enables `BUILD_PYTHON_BINDINGS`, builds position-independent libraries,
and uses shared Boost. Jemalloc is disabled in this configuration to avoid static
TLS failures when importing the extension. Other build directories keep their
normal allocator configuration. Building the native library can take several minutes.

For direct CMake development, use a Python environment with `pybind11` installed:

```sh
cmake -S metagraph -B build-python -DBUILD_PYTHON_BINDINGS=ON
cmake --build build-python --target _pymetagraph_core -j8
cmake --install build-python --component _pymetagraph_core --prefix ./python-install
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
and `extend_chains`. Scoring defaults match the CLI.
With `connect_anchors=False`, a result containing an unresolved graph gap raises
`ValueError`: raw chains cannot use the alignment JSON schema. Keep anchor
connection enabled to obtain inexact alignments. Batch failures drain all workers
before raising, and the loaded graph remains usable.

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
