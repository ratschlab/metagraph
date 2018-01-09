cmake_minimum_required (VERSION 3.6.1)
project(MetaGenGraph)

enable_testing()

string(APPEND CMAKE_CXX_FLAGS " \
  -std=c++14 -Wall -Wextra -Werror \
  -msse4.2 -fopenmp -D_THREAD_SAFE -pthread") # -DDBGDEBUG

# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O2")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O3")

set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)

link_directories(
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
)

include_directories(
  external-libraries/libmaus2/include
  external-libraries/sdsl-lite/include
)

#add_subdirectory()

file(GLOB src_files "*.cpp" "*.h" "*.hpp")

list(FILTER src_files EXCLUDE REGEX ".*\\._.*")
list(FILTER src_files EXCLUDE REGEX ".*\\main.cpp")

add_library(metagraph STATIC ${src_files})
add_executable(metagengraph "main.cpp")
add_dependencies(metagengraph metagraph)

set(METALIBS
  -lhts -lz -lbz2
  -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256
  -lsdsl -ldivsufsort -ldivsufsort64
  -lssl -lcrypto -llzma
)

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -static)
  set(METALIBS ${METALIBS} -lcurl)
endif()

target_link_libraries(
  metagengraph metagraph ${METALIBS}
)


#-------------------
# Unit Tests
#-------------------
# Download and unpack googletest at configure time
if(NOT EXISTS ${CMAKE_BINARY_DIR}/googletest-download)
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(
  ${CMAKE_BINARY_DIR}/googletest-src
  ${CMAKE_BINARY_DIR}/googletest-build
  EXCLUDE_FROM_ALL
)

file(GLOB unit_tests_files "tests/*.cpp")
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories(
  .
)

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files})
add_dependencies(unit_tests metagraph)

target_link_libraries(unit_tests -lgtest_main -lgtest metagraph ${METALIBS})

add_test(NAME unit_tests COMMAND unit_tests)


#-------------------
# Python interfaces
#-------------------
if(PYTHON_INTERFACE)
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  FIND_PACKAGE(PythonLibs)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  SET(CMAKE_SWIG_FLAGS "")

  file(GLOB interface_files "interfaces/*.i")

  SET_SOURCE_FILES_PROPERTIES(${interface_files} PROPERTIES CPLUSPLUS ON)

  SWIG_ADD_LIBRARY(metagraph LANGUAGE python SOURCES ${interface_files})
  SWIG_LINK_LIBRARIES(metagraph metagraph ${METALIBS} ${PYTHON_LIBRARIES})
endif()
