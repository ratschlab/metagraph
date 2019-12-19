include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.6.1)
project(metagraph_utils)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)

findOpenMP()
findFolly()
checkFileSystem() # verify if std::filesystem::temp_directory_path() is available

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_utils STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

# ASIO depends by default on boost/config.hpp, which is not recommended; this works around the issue
target_compile_options(metagraph_utils PRIVATE -DASIO_STANDALONE -Werror -Wall -Wextra)

target_include_directories(metagraph_utils
  PUBLIC
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    .. # this ensures include files in the project are visible relative to the project root
)

find_library(SDSL_LIBRARY sdsl HINTS ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
find_library(LIBMAUS_LIBRARY maus2 HINTS ${EXTERNAL_LIB_DIR}/libmaus2/lib)

target_link_libraries(metagraph_utils
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${FILESYSTEM_LIBRARY}
    ${SDSL_LIBRARY}
    ${LIBMAUS_LIBRARY}
)
