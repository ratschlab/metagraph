# AGENTS.md

## Cursor Cloud specific instructions

### Overview

MetaGraph is a C++ bioinformatics tool for scalable construction of annotated genome graphs and sequence-to-graph alignment. The primary deliverable is a single CLI binary (`metagraph`). It also includes a Python API client and Snakemake workflows (both optional).

### System dependencies (pre-installed)

The build requires: `g++` (13+), `cmake` (3.16+), `make`, `autoconf`, `automake`, `ccache`, and libs: `libboost-all-dev`, `libbz2-dev`, `libdeflate-dev`, `libjemalloc-dev`, `liblzma-dev`, `libzstd-dev`, `python3-venv`.

The default system `c++`/`cc` must point to GCC, not Clang. If they point to Clang, fix with:
```
sudo update-alternatives --set c++ /usr/bin/g++-13
sudo update-alternatives --set cc /usr/bin/gcc-13
```

### Building

Refer to the top-level `Makefile` and `README.md` for build commands. Key steps:

1. `git submodule sync && git submodule update --init --recursive`
2. `mkdir -p metagraph/build && cd metagraph/build`
3. `CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake -DCMAKE_DBG_ALPHABET=DNA -DBUILD_KMC=OFF ..`
4. `make metagraph -j$(nproc)`

The binary is at `metagraph/build/metagraph_DNA` (symlinked as `metagraph/build/metagraph`).

### Lint

There is no separate lint command. The CMake build uses `-Wall -Wextra -Werror`, so a successful build implies lint passes.

### Testing

- **Unit tests**: Build with `make unit_tests -j$(nproc)` in the build dir, run with `./unit_tests`. The full suite takes 5+ minutes; use `--gtest_filter` for subsets.
- **Integration tests**: Run `./integration_tests` in the build dir (requires the test venv created during cmake).
- **Python API tests**: `cd metagraph/api/python && tox` (requires tox installed).

### Gotchas

- `BUILD_KMC=OFF` is the default in the Makefile's `additional_cmake_args` and avoids building the KMC external tool (saves build time). The cmake default is `BUILD_KMC=ON`.
- The cmake step automatically creates a Python venv at `metagraph/build/test_venv` and installs the Python API + test deps into it. This requires `python3-venv` to be installed.
- The full unit test suite is large. For quick verification, filter with `--gtest_filter`.
- Alphabet variants (`DNA`, `DNA5`, `Protein`) produce separate binaries (`metagraph_DNA`, `metagraph_DNA5`, `metagraph_Protein`). The default is DNA.
