cmake_minimum_required (VERSION 3.6.1)
project(MetaGenGraph)

enable_testing()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif()


if(NOT CMAKE_DBG_ALPHABET)
  set(CMAKE_DBG_ALPHABET "DNA")
endif()

if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA_GRAPH")
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA_CASE_SENSITIVE_GRAPH")
else()
  string(APPEND CMAKE_CXX_FLAGS " -D_PROTEIN_GRAPH")
endif()


string(APPEND CMAKE_CXX_FLAGS " \
  -std=c++14 -Wall -Wextra -Werror \
  -msse4.2 -mavx -mavx2 -mfma -mbmi -mbmi2 \
  -fopenmp -D_THREAD_SAFE -pthread") # -DDBGDEBUG -Wconversion -Wsign-conversion


# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O2 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

# Thread sanitizer build type
set(CMAKE_CXX_FLAGS_THREADS "-DNDEBUG -O2 -g -fsanitize=thread")
set(CMAKE_EXE_LINKER_FLAGS_THREADS "-g -fsanitize=thread")
set(CMAKE_SHARED_LINKER_FLAGS_THREADS "-g -fsanitize=thread")

set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)


#-------------------
# KMC k-mer counter
#-------------------
set(KMC_MAIN_DIR "${PROJECT_SOURCE_DIR}/external-libraries/KMC/kmer_counter")
if(NOT DEFINED BUILD_KMC)
set(BUILD_KMC TRUE)
endif()
configure_file(CMakeListsKMC.txt.in KMC/CMakeLists.txt @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/KMC
)
if(result)
  message(FATAL_ERROR "CMake step for KMC failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release -- -j
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/KMC
)
if(result)
  message(FATAL_ERROR "Build step for KMC failed: ${result}")
endif()


link_directories(
  ${CMAKE_BINARY_DIR}/KMC
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
)

include_directories(
  external-libraries/KMC/kmc_api
  external-libraries/libmaus2/include
  external-libraries/sdsl-lite/include
  external-libraries/hopscotch-map/src
  external-libraries/rollinghashcpp
  external-libraries/ips4o
  external-libraries/eigen
  external-libraries/caches/include
  external-libraries/graph_annotation/bloom_annotator
)

add_subdirectory(
  external-libraries/graph_annotation/bloom_annotator
)

set(FOLLY_LIBRARIES FOLLY_LIBRARIES-NOTFOUND)
find_library(FOLLY_LIBRARIES
  NAMES folly
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
  ~/.linuxbrew/lib/
)
if(FOLLY_LIBRARIES)
  string(APPEND CMAKE_CXX_FLAGS " -D_USE_FOLLY")
  set(METALIBS -lfolly)
else()
  message("Folly was not found. Install folly to get most out of the row-based annotator.")
endif()

file(GLOB src_files "*.cpp" "*.h" "*.hpp")

list(FILTER src_files EXCLUDE REGEX ".*\\._.*")
list(FILTER src_files EXCLUDE REGEX ".*main.cpp")

add_library(metagraph STATIC ${src_files})
add_executable(metagengraph "main.cpp")
set_target_properties(metagengraph PROPERTIES OUTPUT_NAME "metagengraph_${CMAKE_DBG_ALPHABET}")

set(METALIBS ${METALIBS}
  -lKMC
  -lhts -lz -lbz2
  -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256
  -lsdsl -ldivsufsort -ldivsufsort64
  -lssl -lcrypto -llzma
  -latomic
)

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -static)
  set(METALIBS ${METALIBS} -lcurl -lssl)
endif()

target_link_libraries(
  metagraph bloom_annotator ${METALIBS}
)

target_link_libraries(
  metagengraph metagraph ${METALIBS}
)

add_custom_target(link_target ALL
  COMMAND ${CMAKE_COMMAND} -E
  create_symlink
    "metagengraph_${CMAKE_DBG_ALPHABET}"
    "metagengraph"
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
