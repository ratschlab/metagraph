cmake_minimum_required(VERSION 3.16)

project(add_subdirectory_example)

set(CMAKE_CXX_STANDARD 17)
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
