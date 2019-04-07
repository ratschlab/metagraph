cmake_minimum_required (VERSION 3.6.1)
project(Metagraph)
include(CheckFunctionExists)
include(CheckLibraryExists)

enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()

set(JSONCPP_WITH_TESTS FALSE)
add_subdirectory(
  external-libraries/jsoncpp
)

if(NOT CMAKE_DBG_ALPHABET)
  set(CMAKE_DBG_ALPHABET "DNA")
endif()

if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA_GRAPH")
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA4")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA4_GRAPH")
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA_CASE_SENSITIVE_GRAPH")
else()
  string(APPEND CMAKE_CXX_FLAGS " -D_PROTEIN_GRAPH")
endif()


string(APPEND CMAKE_CXX_FLAGS " \
  -Wall -Wextra -Werror \
  -msse4.2 \
  -D_THREAD_SAFE -pthread \
  -Wno-unused-parameter \
") # -DDBGDEBUG -Wconversion -Wsign-conversion -Werror=shadow=compatible-local

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " \
    -Wno-exit-time-destructors \
    -Wno-unused-function \
    -Wno-old-style-cast \
    -Wno-atomic-alignment \
    -Wno-unused-parameter \
    -Wno-extra-semi-stmt \
  ")
endif()

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(NOT DEFINED WITH_AVX)
set(WITH_AVX TRUE)
endif()
if(WITH_AVX)
string(APPEND CMAKE_CXX_FLAGS " -mavx -mavx2 -mfma -mbmi -mbmi2")
endif()

if(WITH_MMAP)
message(STATUS "Enabling boost iostreams for mmapped reader")
string(APPEND CMAKE_CXX_FLAGS " -D_USE_MMAP")
set(METALIBS ${METALIBS} -lboost_iostreams)
endif()

if(WITH_DEATH_TESTS)
message(STATUS "Enabling death tests in gtest")
string(APPEND CMAKE_CXX_FLAGS " -D_DEATH_TEST")
endif()

#-------------------
# KMC k-mer counter
#-------------------
set(KMC_MAIN_DIR "${PROJECT_SOURCE_DIR}/external-libraries/KMC/kmer_counter")
if(NOT DEFINED BUILD_KMC)
set(BUILD_KMC TRUE)
endif()
configure_file(CMakeListsKMC.txt.in KMC/CMakeLists.txt @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
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


# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O3 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")

# Link time optimization
if(LINK_OPT)
  string(APPEND CMAKE_CXX_FLAGS " -flto")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -flto")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_AR  "gcc-ar")
    set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_CXX_ARCHIVE_FINISH   true)
  endif()
endif()

# Thread sanitizer build type
set(CMAKE_CXX_FLAGS_THREADS "-DNDEBUG -O2 -g -fsanitize=thread")
set(CMAKE_EXE_LINKER_FLAGS_THREADS "-g -fsanitize=thread")
set(CMAKE_SHARED_LINKER_FLAGS_THREADS "-g -fsanitize=thread")


set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)

link_directories(
  ${CMAKE_BINARY_DIR}/KMC
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_BINARY_DIR}/external-libraries/jsoncpp/src/lib_json
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
)

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/* ${curdir}/*/*)
  SET(dirlist "${curdir}")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")

include_directories(
  external-libraries/KMC/kmc_api
  external-libraries/libmaus2/include
  external-libraries/sdsl-lite/include
  external-libraries/rollinghashcpp
  external-libraries/ips4o
  external-libraries/eigen
  external-libraries/caches/include
  external-libraries/sparsepp
  external-libraries/hopscotch-map/include
  external-libraries/ordered-map/include
  external-libraries/binrel_wt/include
  external-libraries/jsoncpp/include
  external-libraries/asio/asio/include
  external-libraries/cpp_progress_bar
  ${INCLUDE_DIRS}
)

add_library(progress_bar STATIC
  external-libraries/cpp_progress_bar/progress_bar.cpp
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/external-libraries/binrel_wt/cmake/modules")
add_subdirectory(
  external-libraries/binrel_wt
)

set(FOLLY_LIBRARIES FOLLY_LIBRARIES-NOTFOUND)
find_library(FOLLY_LIBRARIES
  NAMES folly
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
    ~/.linuxbrew/lib
)
if(FOLLY_LIBRARIES)
  string(APPEND CMAKE_CXX_FLAGS " -D_USE_FOLLY")
  set(METALIBS ${METALIBS}
    ${FOLLY_LIBRARIES}
    -ldouble-conversion
    -ldl
    -lboost_regex
    -lboost_filesystem
    -lboost_program_options
  )
  if(NOT BUILD_STATIC)
    set(METALIBS ${METALIBS}
      -lglog
      -lgflags
    )
  endif()
  find_library(LIBERTY
    NAMES iberty
    HINTS
    ENV LD_LIBRARY_PATH
    PATHS
      ~/.linuxbrew/lib
  )
  if(LIBERTY)
    set(METALIBS ${METALIBS} ${LIBERTY})
  endif()
else()
  message(WARNING "Folly was not found.\n"
  "Install folly to reduce RAM required by the row-compressed annotator.")
endif()

file(GLOB src_files
  "src/*.cpp"
  "src/*/*.cpp"
  "src/*/*/*.cpp"
)

list(FILTER src_files EXCLUDE REGEX ".*\\._.*")
list(FILTER src_files EXCLUDE REGEX ".*main.cpp")

add_library(metagraph STATIC ${src_files})
add_executable(metagengraph "src/main.cpp")
set_target_properties(metagengraph PROPERTIES OUTPUT_NAME "metagengraph_${CMAKE_DBG_ALPHABET}")

set(METALIBS ${METALIBS}
  -lKMC
  -lhts -lz -lbz2
  -lmaus2
  -lsdsl -ldivsufsort -ldivsufsort64
  -lssl -lcrypto -llzma
  -lbrwt
  -ljsoncpp
  progress_bar
)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(METALIBS ${METALIBS} -latomic)
endif()

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -static)
  set(METALIBS ${METALIBS} -lcurl -lssl)
endif()

# check for std::filesystem::temp_directory_path
include(CheckCXXSourceRuns)
set(CMAKE_REQUIRED_FLAGS " -std=c++17")
set(CMAKE_REQUIRED_LIBRARIES "c++fs")
check_cxx_source_runs("
    #include <iostream>
    #include <filesystem>
    int main() {
        std::cout << std::filesystem::temp_directory_path();
        return 0;
    }
" CPPFS)
unset(CMAKE_REQUIRED_FLAGS)
unset(CMAKE_REQUIRED_LIBRARIES)
if(CPPFS)
  set(METALIBS ${METALIBS} -lc++fs)
else()
  set(CMAKE_REQUIRED_FLAGS " -std=c++17")
  set(CMAKE_REQUIRED_LIBRARIES "stdc++fs")
  check_cxx_source_runs("
      #include <iostream>
      #include <filesystem>
      int main() {
          std::cout << std::filesystem::temp_directory_path();
          return 0;
      }
  " STDCPPFS)
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(STDCPPFS)
    set(METALIBS ${METALIBS} -lstdc++fs)
  else()
    message(FATAL_ERROR "std::filesystem not found")
  endif()
endif()

target_link_libraries(
  metagraph ${METALIBS}
)

target_link_libraries(
  metagengraph metagraph ${METALIBS}
)

add_custom_target(link_target ALL
  COMMAND ${CMAKE_COMMAND} -E
  create_symlink
    "metagengraph_${CMAKE_DBG_ALPHABET}"
    "metagengraph"
)


#-------------------
# Experiments
#-------------------
file(GLOB experiment_files
  "experiments/*.cpp"
)
list(FILTER experiment_files EXCLUDE REGEX ".*\\._.*")

include_directories(
  .
  ${INCLUDE_DIRS}
  experiments
  external-libraries/tclap/include
)

link_directories(
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
  external-libraries/binrel_wt/lib
)
add_executable(run_experiments ${experiment_files})

target_link_libraries(
  run_experiments metagraph ${METALIBS}
)

#-------------------
# install library
#-------------------

file(GLOB header_files
  "src/*.hpp"
  "src/*/*.hpp"
  "src/*/*/*.hpp"
)

install(TARGETS metagraph DESTINATION lib)
install(FILES ${header_files} DESTINATION include/metagraph)

#-------------------
# Unit Tests
#-------------------
if(FALSE) # there is one gtest version attached with sdsl
# Download and unpack googletest at configure time
if(NOT EXISTS ${CMAKE_BINARY_DIR}/googletest-download)
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(
  ${CMAKE_BINARY_DIR}/googletest-src
  ${CMAKE_BINARY_DIR}/googletest-build
  EXCLUDE_FROM_ALL
)
endif()

file(GLOB unit_tests_files
  "tests/*.cpp"
  "tests/*/*.cpp"
  "tests/*/*/*.cpp"
)
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories(
  ${INCLUDE_DIRS}
)

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files})

target_link_libraries(unit_tests gtest_main gtest
  metagraph
  ${METALIBS}
)

add_test(NAME unit_tests COMMAND unit_tests)


#-------------------
# Python interfaces
#-------------------
if(PYTHON_INTERFACE)
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  FIND_PACKAGE(PythonLibs)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  SET(CMAKE_SWIG_FLAGS "")

  file(GLOB interface_files "interfaces/*.i")

  SET_SOURCE_FILES_PROPERTIES(${interface_files} PROPERTIES CPLUSPLUS ON)

  SWIG_ADD_LIBRARY(metagraph LANGUAGE python SOURCES ${interface_files})
  SWIG_LINK_LIBRARIES(metagraph metagraph ${METALIBS} ${PYTHON_LIBRARIES})
endif()
