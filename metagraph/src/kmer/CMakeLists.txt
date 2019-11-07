#Library for extracting, sorting, de-duping and storing kmers for further processing
cmake_minimum_required(VERSION 3.6.1)
project(kmer)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../common common)

file(GLOB src_files "*.cpp" "kmc_counter/*.cpp")

add_library(kmer STATIC ${src_files})
set(EXTERNAL_LIB_DIR "../../external-libraries")

target_include_directories(kmer PUBLIC ${EXTERNAL_LIB_DIR}/sdsl-lite/include)
# this ensures that this libs include files are visible to its deps
target_include_directories(kmer PUBLIC ..)

target_link_libraries(kmer  ${FOLLY_LIBRARIES} -lsdsl)
