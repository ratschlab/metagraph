# pymetagraph

Python bindings for querying an annotated MetaGraph index in-process, without a running
`metagraph server`. They expose the label query behind `metagraph query` for batches of sequences
and return plain `(label, k-mer count)` results. To query a remote server, use the `metagraph`
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

`query` releases the GIL and spreads the sequences over `threads` threads.

## Tests

The tests build a small index with the `metagraph` binary, taken from `$METAGRAPH` or from `PATH`,
and compare the results with `metagraph query`:

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
