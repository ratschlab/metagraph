cmake_minimum_required(VERSION 3.13)
#project(compression_naive_code)
#
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -lc++fs")
#enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIR})
#
#add_executable(compression_naive_code main.cpp)
project(lossless_dbg)
if(NOT DEFINED METALIBS)
    set(METALIBS metagraph sdsl maus2 z hts c++fs)
    string(APPEND CMAKE_CXX_FLAGS " -D_DNA_GRAPH")
endif()

add_library(path_encoding
        utilities.cpp
        path_database_baseline.cpp)

target_link_libraries(path_encoding metagraph ${METALIBS})

add_executable(tests tests.cpp)
add_executable(path_encoder_toolbox main.cpp)


target_link_libraries(tests metagraph ${METALIBS} path_encoding gtest gtest_main)
target_link_libraries(path_encoder_toolbox metagraph ${METALIBS} path_encoding gtest gtest_main )

install(TARGETS path_encoder_toolbox DESTINATION bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(/usr/local/include/metagraph)
#add_subdirectory(playground/compressor)
#add_subdirectory(playground/decompressor)
#add_subdirectory(playground/sampler)
#add_subdirectory(playground/tester)
