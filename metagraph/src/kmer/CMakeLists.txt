#Library for extracting, sorting, de-duping and storing kmers for further processing
cmake_minimum_required(VERSION 3.6.1)
project(kmer)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../common common)
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
message("Executing ${CMAKE_COMMAND} -G \"${CMAKE_GENERATOR}\" ${PROJECT_PARENT_DIR}/external-libraries/KMC/")
execute_process(COMMAND ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR} "${PROJECT_PARENT_DIR}/external-libraries/KMC/"
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

target_link_libraries(kmer  ${FOLLY_LIBRARIES} -lsdsl)
