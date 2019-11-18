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

findOpenMP()
buildProgressBar(EXTERNAL_LIB_DIR progress_bar)
buildMersenneTwister(EXTERNAL_LIB_DIR mersenne_twister)
find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)

enable_testing()
find_package(GTest REQUIRED)

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})

target_include_directories(tests PRIVATE
    ${GTEST_INCLUDE_DIRS}
)
target_link_libraries(tests
    ${GTEST_BOTH_LIBRARIES}
    ${SDSL_LIBRARY}
    kmer
    ${mersenne_twister}
    metagraph_graph_base
    metagraph_succinct
    metagraph_utils
    ${progress_bar}
    z # yes, someone thought that giving zlib a one letter name was a smart idea
    )

