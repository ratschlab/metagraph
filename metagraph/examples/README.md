# How to Use MetaGraph as a C++ Library

Examples demonstrating how to integrate MetaGraph into your C++ project using CMake.

## Integration Methods

Two approaches for integrating MetaGraph:

### 1. [add_subdirectory](subdirectory/)

Simple and direct approach using `add_subdirectory`.

- Best for: Active development, simpler setup
- See: [subdirectory/README.md](subdirectory/README.md)

### 2. [ExternalProject_Add](ExternalProject/)

Automatic download and isolated build.

- Best for: Stricter isolation
- See: [ExternalProject/README.md](ExternalProject/README.md)

## Example Code

[`src/basic_query.cpp`](src/basic_query.cpp) demonstrates:

- Loading a de Bruijn graph and annotations
- Querying sequences against the annotated graph
- Retrieving matching labels

## Test Data

Pre-built graphs and query sequences in [`data/`](data/):

- DNA and Protein graphs with annotations
- Sample query sequences

Both examples include a `make run_example` target for quick testing.
