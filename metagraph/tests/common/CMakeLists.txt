cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../../src/common metagraph_common)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

enable_testing()
findGTest()

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_include_directories(tests
  PRIVATE
    ${GTEST_INCLUDE_DIRS}
    "../../src/"
    ${Boost_INCLUDE_DIRS}
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
)
target_link_libraries(tests
    metagraph_common
    ${GTEST_BOTH_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
    ${SDSL_LIBRARY}
    ${EXTERNAL_LIB_DIR}/libmaus2/lib/libmaus2.a
)

add_executable(file_merger_benchmark file_merger_benchmark.cpp)

add_subdirectory(${PROJECT_PARENT_DIR}/external-libraries/benchmark EXCLUDE_FROM_ALL benchmark)

target_include_directories(file_merger_benchmark
  PRIVATE
    ${EXTERNAL_LIB_DIR}/benchmark/include
    "../../src/"
)

target_link_libraries(file_merger_benchmark
    PRIVATE
    benchmark_main
    benchmark
    metagraph_common
    )

