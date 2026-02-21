cmake_minimum_required(VERSION 3.16)

project(ExternalProject_Add_example)

set(CMAKE_CXX_STANDARD 17)
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
if(DEFINED BUILD_STATIC)
    list(APPEND METAGRAPH_CMAKE_ARGS -DBUILD_STATIC=${BUILD_STATIC})
endif()
if(DEFINED WITH_AVX)
    list(APPEND METAGRAPH_CMAKE_ARGS -DWITH_AVX=${WITH_AVX})
endif()
if(DEFINED WITH_MSSE42)
    list(APPEND METAGRAPH_CMAKE_ARGS -DWITH_MSSE42=${WITH_MSSE42})
endif()

# Download and build MetaGraph as an external project
ExternalProject_Add(metagraph-external
    GIT_REPOSITORY https://github.com/ratschlab/metagraph.git
    GIT_TAG ${METAGRAPH_GIT_TAG}
    SOURCE_SUBDIR metagraph
    CMAKE_ARGS ${METAGRAPH_CMAKE_ARGS}
)

# Get the install directory
ExternalProject_Get_Property(metagraph-external INSTALL_DIR)
set(METAGRAPH_INSTALL_DIR ${INSTALL_DIR})

# Create our example executable
add_executable(external_query ../src/basic_query.cpp)

# On macOS, add Homebrew library paths and detect OpenMP (same as main MetaGraph)
if(APPLE)
    target_link_directories(external_query PRIVATE /usr/local/lib)
    if(EXISTS /opt/homebrew/lib)
        target_link_directories(external_query PRIVATE /opt/homebrew/lib)
    endif()
    target_include_directories(external_query SYSTEM PRIVATE
        /usr/local/include
        /opt/homebrew/include
    )
    if (NOT DEFINED OpenMP_ROOT)
        execute_process(
            COMMAND brew --prefix libomp
            OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(HOMEBREW_LIBOMP_PREFIX)
            set(OpenMP_ROOT "${HOMEBREW_LIBOMP_PREFIX}" CACHE PATH "Path to OpenMP")
        endif()
    endif()
endif()

# Find required system dependencies for MetaGraph
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams)

# Use the installed includes and libraries
target_include_directories(external_query PRIVATE ${METAGRAPH_INSTALL_DIR}/include)
target_link_directories(external_query PRIVATE ${METAGRAPH_INSTALL_DIR}/lib)

# Link MetaGraph and external system dependencies
# metagraph-core.a now bundles all internal dependencies (zlib, sdsl, hts, spdlog, folly, etc.)
# We only need to link against external system dependencies
set(EXTERNAL_LIBS
    metagraph-cli
    metagraph-core
    Boost::iostreams
    OpenMP::OpenMP_CXX
    jemalloc
    deflate
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
