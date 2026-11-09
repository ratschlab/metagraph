# pymetagraph

Python bindings for querying an annotated MetaGraph index in-process, without a running
`metagraph server_query` process. For batches of sequences, they expose the label query behind
`metagraph query`, sequence-to-graph alignment as in `metagraph align`, and the seed-only label
mapping of `metagraph align --align-map-only`, and return plain values. To query a remote server, use the `metagraph`
client package in `metagraph/api/python` instead.

The bindings are experimental.

## Building

Initialize the git submodules (see the top-level README), then build and install from this
directory:

```
pip install .
```

This compiles MetaGraph with `-DBUILD_PYTHON_BINDINGS=ON` through `scikit-build-core`, so it takes
as long as a regular MetaGraph build. CMake options can be passed with the `CMAKE_ARGS`
environment variable, for example `CMAKE_ARGS="-DBoost_DIR=..."` to select a Boost installation.
After a failed build, delete `build/` before retrying.

To rebuild the extension incrementally while working on the C++ side, configure CMake directly
from the repository root, then link the built module into the package:

```
cmake -S metagraph -B build -DBUILD_PYTHON_BINDINGS=ON -DPython_EXECUTABLE=$(which python3)
cmake --build build --target _pymetagraph_core
ln -s $PWD/build/_pymetagraph_core*.so metagraph/api/python/pymetagraph/pymetagraph/
```

and put `metagraph/api/python/pymetagraph` on `PYTHONPATH`.

## Usage

```python
import pymetagraph

graph = pymetagraph.AnnotatedGraph("graph.dbg", "graph.column.annodbg", threads=8)
graph.k        # 31
graph.labels   # ['ref_A', 'ref_B']

hits = graph.query(["ACGT...", "TTGA..."], num_top_labels=1, discovery_fraction=0.1)
# [[('ref_A', 412)], []]   one list of (label, matching k-mers) per sequence
graph.query("ACGT...")     # a single sequence gives a single list
```

`query` has the defaults of `metagraph query`: all labels (`num_top_labels=None`),
`discovery_fraction=0.7` (`--min-kmers-fraction-label`) and `presence_fraction=0.0`.

```python
alignments = graph.align(["ACGT...", "TTGA..."], min_exact_match=0.3)
# [[GraphAlignment(labels=['ref_A'], score=612, cigar='300=', orientation=False,
#                  query_begin=0, query_end=300)], []]   best first, per sequence
graph.align(seqs, labeled=False)   # no labels, as the server's /align endpoint

mappings = graph.map(["ACGT...", "TTGA..."], min_exact_match=0.3)
# [LabelMapping(labels=['ref_A'], num_matches=287, orientation=False,
#               query_begin=0, query_end=300), None]
```

`align` follows the annotation and reports its labels (`metagraph align -a`); with
`labeled=False` it aligns to the whole graph without labels, as the `/align` endpoint of
`metagraph server_query` does. `map` finds and label-filters seeds as labeled alignment does but
does not extend them, and reports the labels whose seeds cover the most of the sequence, which
is enough for accept/reject decisions at a fraction of the cost. Both take `min_exact_match`
(`--align-min-exact-match`, default 0.7), `min_seed_length` (default `min(k, 19)`) and
`max_seed_length`; `align` also takes `max_alternative_alignments` (default 1) and
`max_nodes_per_seq_char` (default 5.0). Unset arguments keep the `metagraph align` defaults.

All three methods release the GIL and spread the sequences over `threads` threads.

## Tests

The tests build a small index with the `metagraph` binary, taken from `$METAGRAPH` or from `PATH`,
and compare the results with `metagraph query` and `metagraph align`:

```
METAGRAPH=/path/to/metagraph pytest tests
```

## Known limitations

- Exactly one annotation file is supported.
- Missing files raise `FileNotFoundError`, but a graph or annotation that MetaGraph cannot load
  ends the Python process, as the loaders exit on errors.
- On Linux, the extension links Boost dynamically and the wheel does not bundle it, so the Boost
  runtime (`libboost_iostreams`) it was built against must be installed where it is imported.
  The wheel is meant for the machine that built it.
- `BUILD_STATIC=ON` cannot be combined with `BUILD_PYTHON_BINDINGS=ON`.
- The extension is built without jemalloc, which can otherwise exhaust the static TLS space when
  the extension is imported.
