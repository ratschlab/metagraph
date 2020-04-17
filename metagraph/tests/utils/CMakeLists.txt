cmake_minimum_required(VERSION 3.6.1)
project(test_utils)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")
include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)

enable_testing()
findGTest()

file(GLOB test_files "*.cpp")

add_library(test_utils ${test_files})

target_include_directories(test_utils
    PRIVATE
    ${GTEST_INCLUDE_DIRS}
)
target_link_libraries(test_utils
    ${GTEST_BOTH_LIBRARIES}
)
