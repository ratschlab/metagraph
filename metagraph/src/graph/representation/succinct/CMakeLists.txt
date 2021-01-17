cmake_minimum_required(VERSION 3.6.1)
project(metagraph_succinct)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

add_subdirectory(${PROJECT_PARENT_DIR}/src/common metagraph_common)
add_subdirectory(${PROJECT_PARENT_DIR}/src/kmer kmer)
add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/representation/base metagraph_graph_base)

findFolly()
checkFileSystem()

file(GLOB src_files "*.cpp" "*.hpp")

add_library(metagraph_succinct STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

target_include_directories(metagraph_succinct
  PUBLIC
    ${PROJECT_PARENT_DIR}/src/
    ${EXTERNAL_LIB_DIR}/rollinchashcpp
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
  PRIVATE
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/cpp_progress_bar
)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)

buildProgressBar(${EXTERNAL_LIB_DIR} progress_bar)
buildMersenneTwister(${EXTERNAL_LIB_DIR} mesenne_twister)

target_link_libraries(metagraph_succinct
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${FILESYSTEM_LIBRARY}
    ${SDSL_LIBRARY}
    mersenne_twister
    progress_bar
    kmer
    metagraph_common
    metagraph_graph_base
)

