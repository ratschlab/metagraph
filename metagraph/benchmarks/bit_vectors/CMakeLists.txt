cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()
enable_testing()
findGTest()

# Check if 128 bit integers are supported
try_compile(MODE_TI
  "${PROJECT_BINARY_DIR}"
  "${EXTERNAL_LIB_DIR}/sdsl-lite/CMakeModules/check_mode_ti.cpp"
)
if(MODE_TI)
  string(APPEND CMAKE_CXX_FLAGS " -DMODE_TI")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(COMMAND sysctl -n machdep.cpu.features
                  OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND sysctl -n machdep.cpu.leaf7_features
                  OUTPUT_VARIABLE _cpuinfo2 OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(_cpuinfo "${_cpuinfo} ${_cpuinfo2}")
endif()

if(WITH_MSSE42 OR (NOT DEFINED WITH_MSSE42 AND _cpuinfo MATCHES "(sse4_2)|(sse4a)|(SSE4.2)"))
string(APPEND CMAKE_CXX_FLAGS " -msse4.2")
else()
message(WARNING "Compiling without msse4.2 instructions!")
endif()

if(WITH_AVX OR (NOT DEFINED WITH_AVX AND _cpuinfo MATCHES "(avx2)|(AVX2)"))
string(APPEND CMAKE_CXX_FLAGS " -mavx -mavx2 -mfma -mbmi -mbmi2")
else()
message(WARNING "Compiling without AVX instructions!")
endif()

message(STATUS "Compile flags used:"
  "${CMAKE_CXX_FLAGS}"
  " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
)

include_directories(
  ${EXTERNAL_LIB_DIR}/benchmark/include
  ${EXTERNAL_LIB_DIR}/DYNAMIC/include
  ${EXTERNAL_LIB_DIR}/sdsl-lite/include
  ${EXTERNAL_LIB_DIR}/hopscotch-map/include
  ${EXTERNAL_LIB_DIR}/ordered-map/include
  ${EXTERNAL_LIB_DIR}/ips4o
  ${EXTERNAL_LIB_DIR}/spdlog/include
  "../../src/"
)
link_directories(
  ${EXTERNAL_LIB_DIR}/sdsl-lite/lib
)

file(GLOB benchmark_files "*.cpp"
  ../../src/common/data_generation.cpp
  ../../src/common/vectors/bit_vector.cpp
  ../../src/common/vectors/bitmap.cpp
  ../../src/common/vectors/bitmap_mergers.cpp
  ../../src/common/vectors/vector_algorithm.cpp
  ../../src/common/serialization.cpp
)

add_subdirectory(${PROJECT_PARENT_DIR}/external-libraries/benchmark EXCLUDE_FROM_ALL benchmark)

add_executable(benchmark_bit_vectors ${benchmark_files})
target_link_libraries(benchmark_bit_vectors
  -lsdsl
  ${OpenMP_CXX_LIBRARIES}
  benchmark_main
  benchmark
)
