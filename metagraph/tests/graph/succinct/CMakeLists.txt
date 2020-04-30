# Build file for the unit-tests in graph/succint
cmake_minimum_required(VERSION 3.6.1)
project(test_succint)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/representation/base metagraph_graph_base)
add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/representation/succinct succinct)
add_subdirectory(${PROJECT_PARENT_DIR}/src/kmer kmer)
add_subdirectory(${PROJECT_PARENT_DIR}/src/common metagraph_common)
add_subdirectory(${PROJECT_PARENT_DIR}/tests/utils test_utils)


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
    ${PROJECT_PARENT_DIR}/
    ${EXTERNAL_LIB_DIR}/googletest/googletest/include}
)
target_link_libraries(tests
    kmer
    metagraph_graph_base
    metagraph_succinct
    metagraph_common
    ${GTEST_BOTH_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
    ${SDSL_LIBRARY}
    ${progress_bar}
    ${mersenne_twister}
    test_utils
    z # yes, someone thought that giving zlib a one letter name was a smart idea
)

