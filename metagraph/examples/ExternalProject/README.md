# Using MetaGraph as a Library With `ExternalProject`

Integrate MetaGraph using `ExternalProject_Add` for automatic download and isolated builds.

## Quick Start

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
make run_example
```

## How It Works

See [`CMakeLists.txt`](CMakeLists.txt) for the complete configuration:

1. Uses `ExternalProject_Add` to download MetaGraph from Git
2. Builds MetaGraph in an isolated build tree
3. Installs headers and libraries to a local prefix
4. Links your executable against installed libraries

## Using in Your Project

Minimal example:

```cmake
include(ExternalProject)

ExternalProject_Add(metagraph-external
    GIT_REPOSITORY https://github.com/ratschlab/metagraph.git
    GIT_TAG master
    SOURCE_SUBDIR metagraph
    CMAKE_ARGS
        -DCMAKE_DBG_ALPHABET=DNA
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Get_Property(metagraph-external INSTALL_DIR)

add_executable(your_app main.cpp)

target_include_directories(your_app PRIVATE ${INSTALL_DIR}/include)
target_link_directories(your_app PRIVATE ${INSTALL_DIR}/lib)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams)

target_link_libraries(your_app PRIVATE metagraph-cli metagraph-core
    Boost::iostreams
    OpenMP::OpenMP_CXX
    jemalloc
    deflate
)

add_dependencies(your_app metagraph-external)
```

**⚠️ CRITICAL for ExternalProject_Add**: Always include `metagraph.hpp` as the **FIRST** header:

```cpp
// metagraph.hpp MUST be included FIRST when using ExternalProject_Add!
// It provides compile-time configuration definitions (_DNA_GRAPH, _USE_FOLLY, etc.)
// that affect how other MetaGraph headers are interpreted.
#include "metagraph.hpp"
#include "seq_io/sequence_io.hpp"
// ... other MetaGraph headers and other external libraryes (e.g., sdsl)
```

See [`CMakeLists.txt`](CMakeLists.txt) for the most up-to-date working version, and more configurations you might need, depending on the system.

## When to Use

**Advantages:**

- Automatic Git download
- Isolated build tree
- Pin to specific versions/commits

**Choose this when** you need stricter build isolation.

For simpler setup, see [add_subdirectory](../subdirectory/README.md).
