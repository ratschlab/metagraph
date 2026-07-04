# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MetaGraph is a C++17 tool for scalable construction of annotated de Bruijn genome graphs and
sequence-to-graph alignment/search, built on succinct data structures. Core engine lives in
`metagraph/src`, exposed through a CLI binary (`metagraph`), a C++ server mode, and a Python
package (`metagraph/workflows`, `metagraph/api/python`).

## Build

The build system is CMake, wrapped by a top-level `Makefile` for convenience (mainly for
docker-based builds). Submodules must be initialized first:

```
git submodule sync && git submodule update --init --recursive
```

`sdsl-lite` additionally needs its own install step before the main build:

```
cd metagraph/external-libraries/sdsl-lite && ./install.sh $PWD
```

Standard local build:

```
mkdir -p metagraph/build && cd metagraph/build
cmake ..
make -j $(($(getconf _NPROCESSORS_ONLN) - 1))       # builds the `metagraph` target
make metagraph -j$(nproc)                            # build just the CLI binary
make unit_tests -j$(nproc)
make integration_tests -j$(nproc)
```

Useful `cmake` configure args (pass as `-D...`):
- `CMAKE_BUILD_TYPE=[Debug|Release|Profile|GProfile]` (default `Release`)
- `CMAKE_DBG_ALPHABET=[DNA|DNA5|DNA_CASE_SENSITIVE|Protein]` (default `DNA`) — selects which
  alphabet the graph/annotation templates are instantiated for; this is a compile-time choice,
  not a runtime flag
- `BUILD_STATIC=[ON|OFF]`, `BUILD_KMC=[ON|OFF]`, `LINK_OPT=[ON|OFF]`, `WITH_AVX=[ON|OFF]`,
  `WITH_MSSE42=[ON|OFF]`

Via the Makefile (mostly for docker workflows): `make build-metagraph env=docker alphabet=Protein`.

Source files are picked up by glob in `metagraph/CMakeLists.txt` (`src/*/*.cpp` etc., several
directories deep) — new `.cpp` files under `src/` or `tests/` are automatically included in the
build without editing CMakeLists.txt; only new subdirectories with their own `CMakeLists.txt`
need to be wired in manually via `add_subdirectory`.

## Tests

Unit tests (GoogleTest), from `metagraph/build`:
```
./unit_tests --gtest_filter="*"
./unit_tests --gtest_filter="DBGSuccinct.*"     # run a subset
```
Test data lives in `metagraph/tests/data` (available to tests via the `TEST_DATA_DIR` macro).

Integration tests (Python, drive the built `metagraph` CLI binary end-to-end), from
`metagraph/build`:
```
./integration_tests --test_filter="*"
./integration_tests --test_filter="test_align*"
./integration_tests --gdb                         # run the binary under gdb on failure
```
Integration test sources are `metagraph/integration_tests/test_*.py`; `main.py` discovers them
by filename glob, so a new file must be named `test_*.py` to be picked up.

There is no separate lint target; `-Wall -Wextra -Werror` is enabled compiler-wide, so warnings
fail the build. Code style follows `metagraph/.clang-format` (LLVM-based, 2-space indent).

## Architecture

### Two build products from one glob
`metagraph-core` (everything under `src/` except `src/cli/**` and `main.cpp`) and
`metagraph-cli` (`src/cli/**`) are separate static libraries linked into the final `metagraph`
executable. `unit_tests` and `benchmarks` link both. Keep this split in mind: code in `cli/`
implements command-line subcommands (`build`, `annotate`, `align`, `query`, `server`, `assemble`,
`clean`, `transform_anno`, `transform_graph`, `stats`, `merge`, `augment` — one `.cpp`/`.hpp` pair
per subcommand in `src/cli/`), while `src/graph`, `src/annotation`, `src/kmer`, `src/common`,
`src/seq_io` are the reusable engine.

### Graph representations (`src/graph/representation/`)
Multiple interchangeable de Bruijn graph backends implement a common `SequenceGraph`/`DeBruijnGraph`
interface (`base/sequence_graph.hpp`):
- `succinct/` — BOSS (`boss.*`) + `dbg_succinct.*`, the primary succinct/scalable representation
- `hash/` — `dbg_hash_ordered`, `dbg_hash_fast`, `dbg_hash_string`, `dbg_sshash` (hash-map-backed)
- `bitmap/` — `dbg_bitmap`, a bitmap-indexed representation
- `masked_graph.*` / `canonical_dbg.*` / `rc_dbg.hpp` — wrappers that mask nodes or add
  canonical/reverse-complement views over any underlying graph

`graph_extensions/` adds optional per-node data (weights, cached first-characters) attachable to
any graph. `annotated_dbg.*` couples a graph with an annotation matrix to answer labeled queries.

### Annotation representations (`src/annotation/`)
`representation/{base,row_compressed,column_compressed,annotation_matrix}` implement the
annotator interface; the actual label-matrix encoding is pluggable via `binary_matrix/`
(`row_flat`, `row_sparse`, `row_disk`, `column_sparse`, `multi_brwt`, `rainbowfish`,
`bin_rel_wt`, `row_diff`, ...). `annotation_converters.*` transforms between these encodings
(this backs the CLI's `transform_anno`, e.g. converting column-compressed to Multi-BRWT via a
linkage file — see README for the two-step clustering + construction workflow).

### Alignment: two implementations coexist mid-migration
`src/graph/alignment/` is the original aligner (`dbg_aligner.*`, `aligner_seeder_methods.*`,
`aligner_extender_methods.*`, `aligner_chainer.*`, `aligner_labeled.*`, `alignment.*`).
`src/graph/alignment_redone/` (`aln_seeder.*`, `aln_match.*`, `aln_cigar.*`, `aln_chainer.hpp`,
`aln_query.hpp`) is a from-scratch rewrite in progress on this branch/work — both are compiled
into the build (CMake globs all of `src/`) and both have corresponding unit test files
(`tests/graph/test_aligner*.cpp` vs `test_aligner_redone.cpp`, `tests/annotation/test_aligner_labeled.cpp`
vs `test_aligner_labeled_redone.cpp`). Check which module a task targets before editing — do not
assume `alignment/` is dead code just because `alignment_redone/` exists, and vice versa.

### Alphabet is a compile-time template parameter
The k-mer/character alphabet (`DNA`, `DNA5`, `DNA_CASE_SENSITIVE`, `Protein`) is selected via the
`CMAKE_DBG_ALPHABET` define (`_DNA_GRAPH` / `_DNA5_GRAPH` / `_DNA_CASE_SENSITIVE_GRAPH` /
`_PROTEIN_GRAPH`), not a runtime option — a given `metagraph` binary only supports one alphabet
(hence separate `metagraph_DNA`, `metagraph_Protein`, etc. binaries/conda packages/docker
variants). `src/kmer/alphabets.hpp` and `kmer_extractor.*` are the places that branch on this.

### External libraries
`metagraph/external-libraries/` vendors ~20 dependencies as git submodules (sdsl-lite, folly,
googletest, htslib, sshash, KMC, Simple-Web-Server, spdlog, DYNAMIC, jsoncpp, etc.) — see
`.gitmodules` for the full list. `sdsl-lite` and `htslib` require their own build steps (see
Build section / `ExternalProject_Add(HTSLIB ...)` in `metagraph/CMakeLists.txt`); the rest build
as part of the main CMake configure.

### Python side
`metagraph/api/python` is the client API package for talking to a running `metagraph server`
instance. `metagraph/workflows` (`metagraph_workflows`) is a separate Snakemake-based pipeline
package for common indexing workflows, with its own `setup.py`/`requirements.txt`.
