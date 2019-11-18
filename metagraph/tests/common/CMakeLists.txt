cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../../src/common metagraph_common)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)
include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()

enable_testing()
find_package(GTest REQUIRED)

set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_compile_options(tests PRIVATE -Werror -Wall -Wextra)
target_include_directories(tests PRIVATE
    ${GTEST_INCLUDE_DIRS}
    "../../src/"
    ${Boost_INCLUDE_DIRS}
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include)
target_link_libraries(tests ${GTEST_BOTH_LIBRARIES} metagraph_common ${SDSL_LIBRARY})

