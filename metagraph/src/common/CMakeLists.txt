include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.6.1)
project(metagraph_common)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)
findOpenMP()
findFolly()
find_package(glog REQUIRED)
checkFileSystem() # check for std::filesystem::temp_directory_path

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_common STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

# ASIO depends by default on boost/config.hpp, which is not recommended; this works around the issue
target_compile_options(metagraph_common PRIVATE -DASIO_STANDALONE -Werror -Wall -Wextra)
target_include_directories(metagraph_common
  PRIVATE
    ${PROJECT_PARENT_DIR}/src/
  PUBLIC
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/nvbio/contrib/htslib
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/rollinghashcpp
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    ${EXTERNAL_LIB_DIR}/spdlog/include
    ${EXTERNAL_LIB_DIR}/zlib
    /usr/local/include  # On Mac the default include directories are not added any longer, so adding manually
)

# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()
set(SPDLOG_BUILD_TESTING FALSE)
add_subdirectory(${EXTERNAL_LIB_DIR}/spdlog spdlog)

add_subdirectory(${EXTERNAL_LIB_DIR}/zlib zlib)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_link_libraries(metagraph_common
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${FILESYSTEM_LIBRARY}
    ${SDSL_LIBRARY}
    glog::glog
    spdlog::spdlog
    z
)
target_compile_definitions(metagraph_common PUBLIC _USE_FOLLY MODE_TI)
message(STATUS "Compile flags used:"
  "${CMAKE_CXX_FLAGS}"
  " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
)
