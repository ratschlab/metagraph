cmake_minimum_required(VERSION 3.6.1)
project(test_succint)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../..)

add_subdirectory(${PROJECT_PARENT_DIR}/src/graph/succinct succinct)
add_subdirectory(${PROJECT_PARENT_DIR}/src/kmer kmer)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

enable_testing()
find_package(GTest REQUIRED)

set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_include_directories(tests PRIVATE
    ${GTEST_INCLUDE_DIRS}
    "../../src/"
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include)
target_link_libraries(tests ${GTEST_BOTH_LIBRARIES} ${SDSL_LIBRARY} kmer metagraph_succinct)

