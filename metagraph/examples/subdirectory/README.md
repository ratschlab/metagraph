# Using MetaGraph as a Library With `add_subdirectory`

Integrate MetaGraph using `add_subdirectory` for simple, direct access to the library.

## Quick Start

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
make run_example
```

## How It Works

See [`CMakeLists.txt`](CMakeLists.txt) for the complete configuration:

1. Add MetaGraph: `add_subdirectory(metagraph/metagraph)`
2. Create executable and link libraries
3. Include paths propagate automatically

## Using in Your Project

Minimal CMakeLists.txt:

```cmake
cmake_minimum_required(VERSION 3.16)
project(YourProject)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(path/to/metagraph/metagraph)

add_executable(your_app main.cpp)
target_link_libraries(your_app PRIVATE metagraph-cli metagraph-core)
```

See [`../src/basic_query.cpp`](../src/basic_query.cpp) for API usage.

## When to Use

**Advantages:**

- Simple setup
- Easy debugging into MetaGraph code
- Fast incremental builds

**Choose this when** you're actively developing with MetaGraph or prefer local source control.

For automatic Git downloads or isolated build tree, see [ExternalProject](../ExternalProject/README.md).
