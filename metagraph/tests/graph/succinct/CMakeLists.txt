# Build file for the unit-tests in graph/succint
cmake_minimum_required(VERSION 3.6.1)
project(test_succint)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/succinct succinct)
add_subdirectory(${PROJECT_PARENT_DIR}/src/kmer kmer)
add_subdirectory(${PROJECT_PARENT_DIR}/src/utils metagraph_utils)
add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/base metagraph_graph_base)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)

enableCCache()
findOpenMP()
buildProgressBar(EXTERNAL_LIB_DIR progress_bar)
buildMersenneTwister(EXTERNAL_LIB_DIR mersenne_twister)
find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)

enable_testing()
findGTest()

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})

target_compile_definitions(tests PUBLIC TEST_DATA_DIR="${PROJECT_PARENT_DIR}/tests/data")
target_include_directories(tests
  PRIVATE
    ${GTEST_INCLUDE_DIRS}
)
target_link_libraries(tests
    kmer
    metagraph_graph_base
    metagraph_succinct
    metagraph_utils
    metagraph_common
    ${GTEST_BOTH_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
    ${SDSL_LIBRARY}
    ${progress_bar}
    ${mersenne_twister}
    z # yes, someone thought that giving zlib a one letter name was a smart idea
)

