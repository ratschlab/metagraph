cmake_minimum_required(VERSION 3.13)
#project(compression_naive_code)
#
set(CMAKE_CXX_STANDARD 17)
set(ENCODER_VERSION 2)
set(PATH_ENCODER_CURRENT path_encoder_toolbox_v${ENCODER_VERSION})
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
        )

target_link_libraries(path_encoding metagraph ${METALIBS})

add_executable(tests tests.cpp)
add_executable(${PATH_ENCODER_CURRENT} main.cpp)

target_link_libraries(tests metagraph ${METALIBS} path_encoding gtest gtest_main)
target_link_libraries(${PATH_ENCODER_CURRENT} metagraph ${METALIBS} path_encoding gtest gtest_main )

macro(install_symlink filepath sympath)
    add_custom_target(link_target_encoder ALL
            COMMAND ${CMAKE_COMMAND} -E
            create_symlink
            ${filepath} ${sympath}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            )
endmacro(install_symlink)

install_symlink(${PATH_ENCODER_CURRENT} path_encoder_toolbox)

install(TARGETS ${PATH_ENCODER_CURRENT} DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/path_encoder_toolbox DESTINATION bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(/usr/local/include/metagraph)
#add_subdirectory(playground/compressor)
#add_subdirectory(playground/decompressor)
#add_subdirectory(playground/sampler)
#add_subdirectory(playground/tester)
