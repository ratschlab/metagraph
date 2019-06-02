cmake_minimum_required(VERSION 3.13)
#project(compression_naive_code)
#
set(CMAKE_CXX_STANDARD 17)
set(ENCODER_VERSION 6)
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
string(TOLOWER ${CMAKE_BUILD_TYPE} BTYPE)
set(PATH_ENCODER_CURRENT path_encoder_toolbox_v${ENCODER_VERSION})
set(PATH_ENCODER_CURRENT_TYPE ${PATH_ENCODER_CURRENT}_${BTYPE})

set(CMAKE_CXX_FLAGS_FASTDEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_EXE_LINKER_FLAGS_FASTDEBUG "-g")
set(CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG "-g")

include_directories(playground/dynamic/include)

## Add profiler into release version
link_libraries(profiler)

set(CMAKE_CXX_FLAGS_GPERF "-O2 -g -lprofiler")
set(CMAKE_EXE_LINKER_FLAGS_GPERF "-g -lprofiler")
set(CMAKE_SHARED_LINKER_FLAGS_GPERF "-g -lprofiler")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -lc++fs")
#enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIR})
#
#add_executable(compression_naive_code main.cpp)
if(NOT DEFINED PROJECT_BINARY_DIR)
    set(ROOT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
else()
    set(ROOT_BINARY_DIR ${PROJECT_BINARY_DIR})
endif()
message("Project binary dir: ${PROJECT_BINARY_DIR}")

project(lossless_dbg)
if(NOT DEFINED METALIBS)
    set(METALIBS metagraph sdsl maus2 z hts c++fs)
    string(APPEND CMAKE_CXX_FLAGS " -D_DNA_GRAPH")
endif()


add_executable(tests tests.cpp)
add_executable(${PATH_ENCODER_CURRENT_TYPE} main.cpp)

target_link_libraries(tests metagraph ${METALIBS} gtest gtest_main)
target_link_libraries(${PATH_ENCODER_CURRENT_TYPE} metagraph ${METALIBS} path_encoding gtest gtest_main )

SET(uid 0)
macro(add_symlink filepath sympath)
    MATH(EXPR uid "${uid}+1")
    add_custom_target(link_target_encoder_${uid} ALL
            COMMAND ${CMAKE_COMMAND} -E
            create_symlink
            ${filepath} ${sympath}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endmacro(add_symlink)

add_symlink(${CMAKE_CURRENT_BINARY_DIR}/${PATH_ENCODER_CURRENT_TYPE} ${ROOT_BINARY_DIR}/path_encoder_toolbox)
add_symlink(${PATH_ENCODER_CURRENT_TYPE} path_encoder_toolbox)
add_symlink(${PATH_ENCODER_CURRENT_TYPE} ${PATH_ENCODER_CURRENT})


install(TARGETS ${PATH_ENCODER_CURRENT_TYPE} DESTINATION bin)
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/path_encoder_toolbox DESTINATION bin)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PATH_ENCODER_CURRENT} DESTINATION bin)
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(/usr/local/include/metagraph)
#add_subdirectory(playground/compressor)
#add_subdirectory(playground/decompressor)
#add_subdirectory(playground/sampler)
#add_subdirectory(playground/tester)
