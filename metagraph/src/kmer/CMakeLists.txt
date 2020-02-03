include_guard(GLOBAL)
#Library for extracting, sorting, de-duping and storing kmers for further processing
cmake_minimum_required(VERSION 3.6.1)
project(kmer)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../common metagraph_common)

SET(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

file(GLOB src_files "*.cpp")

add_library(kmer STATIC ${src_files})

set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

target_include_directories(kmer
  PRIVATE
    ${EXTERNAL_LIB_DIR}/ips4o
  PUBLIC
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    ${EXTERNAL_LIB_DIR}/eigen
    ..
)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_link_libraries(kmer
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${SDSL_LIBRARY}
    metagraph_common
)


if(NOT CMAKE_DBG_ALPHABET)
  set(CMAKE_DBG_ALPHABET "DNA")
endif()

if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
  target_compile_definitions(kmer PUBLIC  -D_DNA_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA5")
  target_compile_definitions(kmer PUBLIC  -D_DNA5_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  target_compile_definitions(kmer PUBLIC  -D_DNA_CASE_SENSITIVE_GRAPH)
else()
  target_compile_definitions(kmer PUBLIC  -D_PROTEIN_GRAPH)
endif()
