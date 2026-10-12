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
findFolly()
find_package(glog REQUIRED)
checkFileSystem() # check for std::filesystem::temp_directory_path

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_common STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")
buildProgressBar(${EXTERNAL_LIB_DIR} PROGRESS_BAR_LIBRARY)

# ASIO depends by default on boost/config.hpp, which is not recommended; this works around the issue
target_compile_options(metagraph_common PRIVATE -DASIO_STANDALONE -Werror -Wall -Wextra)
target_include_directories(metagraph_common
  PUBLIC
    ${PROJECT_PARENT_DIR}/src/
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/rollinghashcpp
    /usr/local/include  # On Mac the default include directories are not added any longer, so adding manually
)
target_include_directories(metagraph_common SYSTEM PUBLIC
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    ${EXTERNAL_LIB_DIR}/simde-no-tests
)

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
    ${PROGRESS_BAR_LIBRARY}
    ${FILESYSTEM_LIBRARY}
    OpenMP::OpenMP_CXX
    glog::glog
    spdlog::spdlog
    ZLIB::ZLIB
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
