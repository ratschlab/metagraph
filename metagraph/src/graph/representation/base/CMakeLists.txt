# Build file for the basic graph classes
include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.6.1)
project(metagraph_graph_base)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)

add_subdirectory(${PROJECT_PARENT_DIR}/src/common metagraph_common)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
buildProgressBar(${EXTERNAL_LIB_DIR} progress_bar)


file(GLOB src_files "*.cpp" "*.hpp")

add_library(metagraph_graph_base STATIC ${src_files})

# the interface of this library depends on stdlib only, so all includes can be private!
target_include_directories(metagraph_graph_base
  PRIVATE
    ${PROJECT_PARENT_DIR}/src/
    ${EXTERNAL_LIB_DIR}/cpp_progress_bar
    ${EXTERNAL_LIB_DIR}/rollinchashcpp
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
)

target_link_libraries(metagraph_graph_base
  PUBLIC
    ${SDSL_LIBRARY}
    ${mersenne_twister}
    progress_bar
    metagraph_common
)

