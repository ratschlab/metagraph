include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.19)
project(metagraph_common)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()
find_package(Boost REQUIRED)
find_package(ZLIB REQUIRED)
buildLogging()
buildCommonDependencies()
findFolly()
find_package(glog REQUIRED)
checkFileSystem() # check for std::filesystem::temp_directory_path

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_common STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")
buildProgressBar(${EXTERNAL_LIB_DIR} PROGRESS_BAR_LIBRARY)
buildMersenneTwister(${EXTERNAL_LIB_DIR} MERSENNE_TWISTER_LIBRARY)

target_compile_options(metagraph_common PRIVATE -Werror -Wall -Wextra)
target_include_directories(metagraph_common PUBLIC ${PROJECT_PARENT_DIR}/src)

# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()

target_link_libraries(metagraph_common
  PUBLIC
    folly
    metagraph-common-dependencies
    ${PROGRESS_BAR_LIBRARY}
    ${MERSENNE_TWISTER_LIBRARY}
    ${FILESYSTEM_LIBRARY}
    OpenMP::OpenMP_CXX
    glog::glog
    spdlog::spdlog
    ZLIB::ZLIB
  PRIVATE
    Boost::headers
)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
  target_link_libraries(metagraph_common PUBLIC atomic)
endif()
target_compile_definitions(metagraph_common PUBLIC _USE_FOLLY MODE_TI)
set(_USE_FOLLY 1)
set(MODE_TI 1)
set(SPDLOG_FMT_EXTERNAL 1)
configure_file(compile_definitions.hpp.in
    ${PROJECT_BINARY_DIR}/gen/include/common/compile_definitions.hpp @ONLY)
target_include_directories(metagraph_common PUBLIC ${PROJECT_BINARY_DIR}/gen/include)
message(STATUS "Compile flags used:"
  "${CMAKE_CXX_FLAGS}"
  " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
)
