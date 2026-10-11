cmake_minimum_required(VERSION 3.16)

project(ExternalProject_Add_example)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set MetaGraph version
set(METAGRAPH_GIT_TAG "master" CACHE STRING "Git tag/branch/commit for MetaGraph")

# Set alphabet (DNA or Protein)
if(NOT DEFINED CMAKE_DBG_ALPHABET)
    set(CMAKE_DBG_ALPHABET "DNA" CACHE STRING "Alphabet for MetaGraph (DNA or Protein)")
endif()

# Include ExternalProject module
include(ExternalProject)

# Forward optional build flags to MetaGraph if they're set
set(METAGRAPH_CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_DBG_ALPHABET=${CMAKE_DBG_ALPHABET}
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DMETAGRAPH_BUILD_EXECUTABLE=OFF
)

# Forward optional flags if they're defined
foreach(variable BUILD_STATIC WITH_AVX WITH_MSSE42 Boost_USE_STATIC_LIBS
                 CMAKE_PREFIX_PATH OpenMP_ROOT ICU_ROOT)
    if(DEFINED ${variable})
        string(REPLACE ";" "|" value "${${variable}}")
        list(APPEND METAGRAPH_CMAKE_ARGS "-D${variable}=${value}")
    endif()
endforeach()

# Download and build MetaGraph as an external project
ExternalProject_Add(metagraph-external
    GIT_REPOSITORY https://github.com/ratschlab/metagraph.git
    GIT_TAG ${METAGRAPH_GIT_TAG}
    SOURCE_SUBDIR metagraph
    LIST_SEPARATOR |
    CMAKE_ARGS ${METAGRAPH_CMAKE_ARGS}
)

# Get the install directory
ExternalProject_Get_Property(metagraph-external INSTALL_DIR)
set(METAGRAPH_INSTALL_DIR ${INSTALL_DIR})

# Create our example executable
add_executable(external_query ../src/basic_query.cpp)

# Find required system dependencies for MetaGraph
include(${CMAKE_CURRENT_LIST_DIR}/../../CMakeListsHelpers.txt)
findOpenMP()
find_package(Boost REQUIRED COMPONENTS iostreams)
find_library(JEMALLOC_LIBRARY NAMES jemalloc)
find_library(LIBDEFLATE_LIBRARY NAMES deflate)

# Use the installed includes and libraries
target_include_directories(external_query PRIVATE ${METAGRAPH_INSTALL_DIR}/include)
target_link_directories(external_query PRIVATE ${METAGRAPH_INSTALL_DIR}/lib)

# Link MetaGraph and external system dependencies
# metagraph-core.a now bundles all internal dependencies (zlib, hts, spdlog, folly, etc.)
# We only need to link against external system dependencies
set(EXTERNAL_LIBS
    metagraph-cli
    metagraph-core
    Boost::iostreams
    OpenMP::OpenMP_CXX
    ${JEMALLOC_LIBRARY}
    ${LIBDEFLATE_LIBRARY}
)

# Add atomic library only for GNU compilers or non-macOS systems
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    list(APPEND EXTERNAL_LIBS atomic)
endif()

target_link_libraries(external_query PRIVATE ${EXTERNAL_LIBS})

# Ensure MetaGraph is built before our executable
add_dependencies(external_query metagraph-external)

# Add a custom target to run the example with test data
add_custom_target(run_example
    COMMAND external_query
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/graphs/test_${CMAKE_DBG_ALPHABET}_graph.dbg
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/graphs/test_${CMAKE_DBG_ALPHABET}_graph.column.annodbg
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/test_${CMAKE_DBG_ALPHABET}_query.fa
    DEPENDS external_query
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running external_query example with ${CMAKE_DBG_ALPHABET} test data"
)
