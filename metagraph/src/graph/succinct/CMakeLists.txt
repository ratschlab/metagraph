cmake_minimum_required(VERSION 3.6.1)
project(metagraph_succinct)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE 1)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

add_subdirectory(${PROJECT_PARENT_DIR}/src/common metagraph_common)
add_subdirectory(${PROJECT_PARENT_DIR}/src/utils metagraph_utils)
add_subdirectory(${PROJECT_PARENT_DIR}/src/kmer kmer)

findFolly()
checkFileSystem()

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_succinct STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

# ASIO depends by default on boost/config.hpp, which is not recommended; this works around the issue
target_compile_options(metagraph_succinct PRIVATE -DASIO_STANDALONE -Werror -Wall -Wextra)

target_include_directories(metagraph_succinct
  PUBLIC
    ${PROJECT_PARENT_DIR}/src/
    ${Boost_INCLUDE_DIRS}
    ${EXTERNAL_LIB_DIR}/cpp_progress_bar
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    .. # this ensures that this libs include files are visible to its deps as "#include common/*"
)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)

target_link_libraries(metagraph_succinct
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${FILESYSTEM_LIBRARY}
    ${SDSL_LIBRARY}
    kmer
    metagraph_common
    metagraph_utils
)
