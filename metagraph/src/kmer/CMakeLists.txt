#Library for extracting, sorting, de-duping and storing kmers for further processing
cmake_minimum_required(VERSION 3.6.1)
project(kmer)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../common metagraph_common)
SET(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)

#-------------------
# KMC k-mer counter
#-------------------
set(KMC_MAIN_DIR "${PROJECT_PARENT_DIR}/external-libraries/KMC/kmer_counter")
if(NOT DEFINED BUILD_KMC)
  set(BUILD_KMC TRUE)
endif()
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/KMC)
configure_file(${PROJECT_PARENT_DIR}/CMakeListsKMC.txt.in ${PROJECT_PARENT_DIR}/external-libraries/KMC/CMakeLists.txt
        @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR} "${PROJECT_PARENT_DIR}/external-libraries/KMC/"
        -DCMAKE_CXX_FLAGS="-w"
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
)
#-------------------
# End KMC k-mer counter
#-------------------

# unset OpenMP_CXX_FLAGS. If OpenMP loading fails once, it will be set to NOT_FOUND
# and happily cached, thus messing up the compilation flags forever
unset(OpenMP_CXX_FLAGS CACHE)
if ((CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
  find_package(OPENMP REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


file(GLOB src_files "*.cpp" "kmc_counter/*.cpp")

add_library(kmer STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

target_include_directories(kmer
        PRIVATE
        ${EXTERNAL_LIB_DIR}/KMC/kmc_api
        ${EXTERNAL_LIB_DIR}/ips4o
        PUBLIC
        ${EXTERNAL_LIB_DIR}/sdsl-lite/include
        ..) # this ensures that this libs include files are visible to its deps

target_link_libraries(kmer  ${FOLLY_LIBRARIES} -lsdsl metagraph_common)


if(NOT CMAKE_DBG_ALPHABET)
  set(CMAKE_DBG_ALPHABET "DNA")
endif()

if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
  target_compile_options(kmer PRIVATE  -D_DNA_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA5")
  target_compile_options(kmer PRIVATE  -D_DNA5_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  target_compile_options(kmer PRIVATE  -D_DNA_CASE_SENSITIVE_GRAPH)
else()
  target_compile_options(kmer PRIVATE  -D_PROTEIN_GRAPH)
endif()

