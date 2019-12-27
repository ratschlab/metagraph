include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.6.1)
project(metagraph_utils)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PARENT_DIR ${PROJECT_SOURCE_DIR}/../../..)

include(${PROJECT_PARENT_DIR}/CMakeListsHelpers.txt)

checkFileSystem() # verify if std::filesystem::temp_directory_path() is available

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_utils STATIC ${src_files})
set(EXTERNAL_LIB_DIR "${PROJECT_PARENT_DIR}/external-libraries")

target_compile_options(metagraph_utils PRIVATE -Werror -Wall -Wextra)

target_link_libraries(metagraph_utils
  PUBLIC
    ${FILESYSTEM_LIBRARY}
)
