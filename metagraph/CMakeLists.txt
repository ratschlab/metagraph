cmake_minimum_required (VERSION 2.6.2)
project(MetaGenGraph)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
list(APPEND CMAKE_CXX_FLAGS "\
  -std=c++11 -DNDEBUG -Wall -O3 \
  -msse4.2 -fopenmp -D_THREAD_SAFE -pthread \
  -W -Wno-long-long -Wno-variadic-macros")

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
add_executable(metagengraph ${src_files})

set(CMAKE_EXE_LINKER_FLAGS "-static")
target_link_libraries(
  metagengraph
  -lhts -lz -lbz2
  -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256
  -lsdsl -ldivsufsort -ldivsufsort64
  -lcurl -lssl -lcrypto -llzma
)
