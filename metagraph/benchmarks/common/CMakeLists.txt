cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()
enable_testing()
findGTest()


file(GLOB test_files "*.cpp")

add_subdirectory(${PROJECT_PARENT_DIR}/external-libraries/benchmark EXCLUDE_FROM_ALL benchmark)
add_subdirectory(${PROJECT_PARENT_DIR}/src/common metagraph_common)

add_executable(benchmark_file_merger benchmark_file_merger.cpp)
target_include_directories(benchmark_file_merger
  PRIVATE
    ${EXTERNAL_LIB_DIR}/benchmark/include
    "../../src/"
)
target_link_libraries(benchmark_file_merger
  PRIVATE
    benchmark_main
    benchmark
    metagraph_common
)

add_executable(benchmark_chunked_wait_queue benchmark_chunked_wait_queue.cpp)
target_include_directories(benchmark_chunked_wait_queue
  PRIVATE
    ${EXTERNAL_LIB_DIR}/benchmark/include
    "../../src/"
)
target_link_libraries(benchmark_chunked_wait_queue
  PRIVATE
    benchmark_main
    benchmark
    metagraph_common
)

add_executable(benchmark_elias_fano benchmark_elias_fano.cpp)
target_include_directories(benchmark_elias_fano
    PRIVATE
    ${EXTERNAL_LIB_DIR}/benchmark/include
    "../../src/"
    )
target_link_libraries(benchmark_elias_fano
    PRIVATE
    benchmark_main
    benchmark
    metagraph_common
    )
