cmake_minimum_required(VERSION 3.19)

project(add_subdirectory_example)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add metagraph as a subdirectory (from git submodule)
add_subdirectory(metagraph/metagraph)

# Create example executables
# Create executable from example code
add_executable(basic_query ../src/basic_query.cpp)

# Link against metagraph libraries - includes are inherited automatically
target_link_libraries(basic_query PRIVATE
    metagraph-cli
    metagraph-core
)

# Coordinate context example (requires DNA, coordinates, and a .seqs mapping).
add_executable(context_cursor ../src/context_cursor.cpp)
target_link_libraries(context_cursor PRIVATE metagraph-cli metagraph-core)

# Set alphabet (DNA or Protein)
if(NOT DEFINED CMAKE_DBG_ALPHABET)
    set(CMAKE_DBG_ALPHABET "DNA" CACHE STRING "Alphabet for MetaGraph (DNA or Protein)")
endif()

# Add a custom target to run the example with test data
add_custom_target(run_example
    COMMAND basic_query
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/graphs/test_${CMAKE_DBG_ALPHABET}_graph.dbg
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/graphs/test_${CMAKE_DBG_ALPHABET}_graph.column.annodbg
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/test_${CMAKE_DBG_ALPHABET}_query.fa
    DEPENDS basic_query
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running basic_query example with ${CMAKE_DBG_ALPHABET} test data"
)

# Reuse the CLI entry point to build fixtures, without enabling MetaGraph's
# full unit/integration/benchmark setup in this library-consumer example.
if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
    add_executable(context_fixture_cli EXCLUDE_FROM_ALL metagraph/metagraph/src/main.cpp)
    target_link_libraries(context_fixture_cli PRIVATE metagraph-cli metagraph-core)
    add_custom_target(run_context_example
        COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test_context_cursor.py
            $<TARGET_FILE:context_cursor> $<TARGET_FILE:context_fixture_cli>
        DEPENDS context_cursor context_fixture_cli
        COMMENT "Checking context requests on fresh heap/mmap coordinate indexes"
    )
    add_dependencies(run_example run_context_example)
endif()
