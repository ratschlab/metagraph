# pymetagraph

A [Readfish](https://github.com/LooseLab/readfish) `Aligner` plugin backed by a prebuilt
metagraph graph + annotation index. By default (`method="query"`) it reports label-presence hits
(which reference sequence(s) a read's k-mers match, above a discovery/presence-fraction threshold)
rather than real base-level alignments, so `r_st`/`r_en`/`strand` on the returned alignment
objects are placeholder values. Set `method="align"` for real seed-and-extend alignment with
actual coordinates and orientation instead — see [Usage](#usage) below.

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
require every user of this plugin to remember `LD_PRELOAD=libjemalloc.so.2` before launching
Readfish. Since this plugin is aimed at non-power-users, it trades away folly's jemalloc-backed
allocator optimizations for an extension that just imports cleanly with no extra steps. This only
affects `_pymetagraph_core` — the regular static `metagraph` CLI binary still picks up jemalloc
normally when built without `BUILD_PYTHON_BINDINGS`.

## Usage

Two unrelated query methods are available via the `method` setting:

- `"query"` (default): label-presence hits via `AnnotatedDBG::get_top_labels` — the same
  primitive backing the CLI's `query` command. Fast and tolerant of sequencing errors, but
  `r_st`/`r_en`/`strand` are dummy values (`0`, `len(seq)`, `+1`); only `ctg` (the matched label)
  is meaningful.
  ```toml
  [mapper_settings.pymetagraph]
  input = "/path/to/graph.dbg"
  annotator = "/path/to/graph.column.annodbg"
  threads = 4
  method = "query"
  num_top_labels = 1
  discovery_fraction = .4
  presence_fraction = .1
  ```
- `"align"`: real seed-and-extend alignment via metagraph's alignment_redone module — the same
  code path backing the `metagraph align` CLI command and the server's `/align` endpoint (what
  the old HTTP client's `GraphClient.align` calls). Slower, but `r_st`/`r_en`/`strand` reflect the
  actual matched region and orientation of the read.
  ```toml
  [mapper_settings.pymetagraph]
  input = "/path/to/graph.dbg"
  annotator = "/path/to/graph.column.annodbg"
  threads = 4
  method = "align"
  seed_length = 19
  max_alternative_alignments = 1
  max_num_nodes_per_seq_char = 10.0
  min_exact_match = .7
  connect_anchors = true
  extend_chains = true
  ```

## Known limitations

- `r_st`/`r_en`/`strand` are always dummy values (`0`, `len(seq)`, `+1`) under `method="query"`;
  only `ctg` (the matched label) is meaningful there. `method="align"` reports real values.
- Exactly one annotator file is supported.
- Malformed configuration may, in rare cases, abort the process rather than raise a Python
  exception, since metagraph's underlying CLI config parser calls `exit()` on some invalid inputs.
  `Aligner.validate()` checks the common cases (missing files) up front.
