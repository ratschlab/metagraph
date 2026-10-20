# pymetagraph

Python bindings for querying an annotated metagraph index in-process, without a running
`metagraph server`. They expose the label query behind `metagraph query` (and the server's
`/search` endpoint) for batches of sequences and return plain `(label, k-mer count)` results, so
callers decide how to turn them into classifications. For talking to a remote server, use the
`metagraph` client package in `metagraph/api/python` instead.

## Building

Unlike a typical small pybind11 project, this is not an instant pip install: metagraph needs its
git submodules checked out and takes a while to compile from scratch. But it *is* a plain
`pip install`-driven build — pip's normal build isolation fetches `scikit-build-core` and
`pybind11` (this project's `[build-system].requires`) into a throwaway env on its own, so you
don't need to pre-install anything Python-side yourself. Recommended workflow:

1. From the repo root, make sure submodules are initialized (see the top-level README).
2. `sdsl-lite` is pre-built separately (via its own `install.sh`, not by the main CMake project)
   and the resulting `external-libraries/sdsl-lite/lib/libsdsl.a` is **not** position-independent
   by default, which breaks the final link into a `.so` even with `BUILD_PYTHON_BINDINGS=ON`. It
   only needs to be rebuilt once (PIC static objects link fine into the normal static `metagraph`
   binary too, so this doesn't need to be undone for the default build):
   ```
   cd metagraph/external-libraries/sdsl-lite/build
   ./clean.sh
   cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
         -DCMAKE_INSTALL_PREFIX=$(cd .. && pwd) ..
   make -j$(nproc) && make install
   ```
3. Build and install — this is the only step most users need, and it never requires running
   `cmake` yourself:
   ```
   cd metagraph/api/python/pymetagraph
   pip install -e .
   ```
   Under the hood, pip's `scikit-build-core` build backend invokes CMake for you (with
   `-DBUILD_PYTHON_BINDINGS=ON` already wired in via this project's `pyproject.toml`) inside a
   throwaway build environment, so you don't need CMake installed in your active environment or
   any familiarity with its flags.

   If you're in a conda environment that ships its own (shared-only) Boost, CMake may pick that up
   ahead of the system Boost metagraph needs (static libs), failing with `boost_iostreams_FOUND to
   FALSE`. Since pip's build backend is the one invoking CMake here, pass the override through the
   standard `CMAKE_ARGS` environment variable, which `scikit-build-core` forwards to it:
   ```
   CMAKE_ARGS="-DBoost_DIR=/usr/lib/x86_64-linux-gnu/cmake/Boost-1.83.0" pip install -e .
   ```
   (adjust the path/version to whatever `find /usr/lib -name 'BoostConfig.cmake'` reports on your
   system).

### Faster iteration on the C++ side (advanced, driving CMake directly)

The steps above are a *different, independent path* from what follows: `pip install -e .` never
invokes `cmake` on your command line, it happens inside pip's own build backend. The workflow
below is only for developers iterating on the C++ binding itself, who want incremental rebuilds
without re-running pip's build backend (and its temp-env setup) on every change. It talks to CMake
directly, so **it needs `-D` flags spelled out explicitly** — none of `pyproject.toml`'s settings
or pip's `CMAKE_ARGS` env var apply here, since pip isn't involved:
```
cmake -S metagraph -B metagraph/api/python/pymetagraph/build -DBUILD_PYTHON_BINDINGS=ON
# add -DBoost_DIR=... too if you hit the same conda/Boost conflict as above
cmake --build metagraph/api/python/pymetagraph/build --target _pymetagraph_core -j$(nproc)
cd metagraph/api/python/pymetagraph
pip install -e . --no-build-isolation   # requires scikit-build-core + pybind11 pre-installed
```

### Note: jemalloc is disabled for this extension

`BUILD_PYTHON_BINDINGS=ON` also sets `CMAKE_DISABLE_FIND_PACKAGE_Jemalloc=ON`, so
`_pymetagraph_core` is built without jemalloc even if it's installed on the system. This is
deliberate: dlopen()-ing a library linked against jemalloc can fail at import time with `cannot
allocate memory in static TLS block` (a glibc static-TLS-space limitation), which would otherwise
require every user to remember `LD_PRELOAD=libjemalloc.so.2` before starting Python. It trades
away folly's jemalloc-backed allocator optimizations for an extension that just imports cleanly
with no extra steps. This only
affects `_pymetagraph_core` — the regular static `metagraph` CLI binary still picks up jemalloc
normally when built without `BUILD_PYTHON_BINDINGS`.

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

## Known limitations

- Exactly one annotation file is supported.
- Graphs and annotations must be readable by the metagraph version the bindings are built from.
- The extension is built without jemalloc (see above).
- `pip install` reuses `build/`; after a failed build, delete it before retrying, since its CMake
  cache can point to pip's deleted build environment.
