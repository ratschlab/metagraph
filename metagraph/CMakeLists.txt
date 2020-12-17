cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0015 NEW)

project(MetaGraph)

include(CheckFunctionExists)
include(CheckLibraryExists)
include(CMakeListsHelpers.txt)

enable_testing()
include(GoogleTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_VERBOSE_MAKEFILE 1)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

function(message)
  list(GET ARGV 0 MessageType)
  if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
    list(REMOVE_AT ARGV 0)
    _message(${MessageType} "${BoldRed}" ${ARGV} "${ColourReset}")
  elseif(MessageType STREQUAL WARNING)
    list(REMOVE_AT ARGV 0)
    _message(${MessageType} "${BoldYellow}" ${ARGV} "${ColourReset}")
  elseif(MessageType STREQUAL AUTHOR_WARNING)
    list(REMOVE_AT ARGV 0)
    _message(${MessageType} "${BoldCyan}" ${ARGV} "${ColourReset}")
  else()
    _message(${ARGV})
  endif()
endfunction()

enableCCache()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

if(${CMAKE_BUILD_TYPE} MATCHES RELEASE)
  message(STATUS "Compile in ${Green}Release${ColourReset} mode")
  set(CMAKE_VERBOSE_MAKEFILE 0)
elseif(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  message(STATUS "Compile in ${Red}Debug${ColourReset} mode")
else()
  message(STATUS "Compile in ${Blue}${CMAKE_BUILD_TYPE}${ColourReset} mode")
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
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA5")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA5_GRAPH")
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  string(APPEND CMAKE_CXX_FLAGS " -D_DNA_CASE_SENSITIVE_GRAPH")
else()
  string(APPEND CMAKE_CXX_FLAGS " -D_PROTEIN_GRAPH")
endif()


string(APPEND CMAKE_CXX_FLAGS
  " -Wall -Wextra -Werror"
  " -D_THREAD_SAFE -pthread"
  # -DDBGDEBUG -Wconversion -Wsign-conversion -Werror=shadow=compatible-local
)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " \
    -Wno-exit-time-destructors \
    -Wno-unused-function \
    -Wno-old-style-cast \
    -Wno-atomic-alignment \
    -Wno-extra-semi-stmt \
  ")
endif()

if(APPLE)
  # add path to system libraries but avoid cross-compilation errors
  include_directories(SYSTEM /usr/local/include)
  link_directories(/usr/local/lib)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-Wno-poison-system-directories SUPPORT_NO_POISON_SYS_DIR)
  if(SUPPORT_NO_POISON_SYS_DIR)
    add_compile_options(-Wno-poison-system-directories)
  endif()
endif()

# unset OpenMP_CXX_FLAGS. If OpenMP loading fails once, it will be set to NOT_FOUND
# and happily cached, thus messing up the compilation flags forever
unset(OpenMP_CXX_FLAGS CACHE)
if (CMAKE_CXX_COMPILER_ID MATCHES AppleClang)
  find_package(OPENMP REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Check if 128 bit integers are supported
try_compile(MODE_TI
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/external-libraries/sdsl-lite/CMakeModules/check_mode_ti.cpp"
)
if(MODE_TI)
  string(APPEND CMAKE_CXX_FLAGS " -DMODE_TI -mcx16")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(COMMAND sysctl -n machdep.cpu.features
                  OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND sysctl -n machdep.cpu.leaf7_features
                  OUTPUT_VARIABLE _cpuinfo2 OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(_cpuinfo "${_cpuinfo} ${_cpuinfo2}")
endif()

if(WITH_MSSE42 OR (NOT DEFINED WITH_MSSE42 AND _cpuinfo MATCHES "(sse4_2)|(sse4a)|(SSE4.2)"))
string(APPEND CMAKE_CXX_FLAGS " -msse4.2")
else()
message(WARNING "Compiling without msse4.2 instructions!")
endif()

if(WITH_AVX OR (NOT DEFINED WITH_AVX AND _cpuinfo MATCHES "(avx2)|(AVX2)"))
string(APPEND CMAKE_CXX_FLAGS " -mavx -mavx2 -mfma -mbmi -mbmi2")
else()
message(WARNING "Compiling without AVX instructions!")
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


# Profile with gprof
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O2 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Profile with gperftools
set(CMAKE_CXX_FLAGS_GPROFILE "-fno-omit-frame-pointer -DNDEBUG -O2 -g")
set(CMAKE_EXE_LINKER_FLAGS_GPROFILE "-g")
set(CMAKE_SHARED_LINKER_FLAGS_GPROFILE "-g")

set(PROFILER_LIBRARIES PROFILER_LIBRARIES-NOTFOUND)
find_library(PROFILER_LIBRARIES
  NAMES profiler
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
    ~/.linuxbrew/lib
)
if(${CMAKE_BUILD_TYPE} MATCHES GPROFILE)
  if(PROFILER_LIBRARIES)
    set(METALIBS ${METALIBS} -lprofiler)
  else()
    message(WARNING "Could not find -lprofiler installed")
  endif()
endif()

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
set(CMAKE_CXX_FLAGS_THREADS "-O2 -g -fsanitize=thread")
set(CMAKE_EXE_LINKER_FLAGS_THREADS "-g -fsanitize=thread")
set(CMAKE_SHARED_LINKER_FLAGS_THREADS "-g -fsanitize=thread")

link_directories(
  ${CMAKE_BINARY_DIR}/KMC
  ${CMAKE_BINARY_DIR}/zlib
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_BINARY_DIR}/external-libraries/jsoncpp/src/lib_json
  external-libraries/sdsl-lite/lib
)

include_directories(
  external-libraries/KMC/kmc_api
  external-libraries/DYNAMIC/include
  external-libraries/sdsl-lite/include
  external-libraries/ips4o
  external-libraries/eigen
  external-libraries/caches/include
  external-libraries/hopscotch-map/include
  external-libraries/ordered-map/include
  external-libraries/binrel_wt/include
  external-libraries/jsoncpp/include
  external-libraries/asio/asio/include
  external-libraries/cpp_progress_bar
  external-libraries/nvbio/contrib/priority-deque
  external-libraries/rollinghashcpp
  external-libraries/spdlog/include
  external-libraries/zlib
  external-libraries/Simple-Web-Server
  external-libraries/folly
  ${PROJECT_SOURCE_DIR}/src
)

# use external fmt library in spdlog
set(SPDLOG_FMT_EXTERNAL TRUE)
set(SPDLOG_BUILD_TESTING FALSE)

# makes ASIO not depend on boost
string(APPEND CMAKE_CXX_FLAGS " -DASIO_STANDALONE -DUSE_STANDALONE_ASIO")

add_subdirectory(external-libraries/fmt EXCLUDE_FROM_ALL)

# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()
add_subdirectory(external-libraries/spdlog)
add_subdirectory(external-libraries/zlib zlib)

add_library(progress_bar STATIC
  external-libraries/cpp_progress_bar/progress_bar.cpp
)

add_library(mersenne_twister STATIC
  external-libraries/rollinghashcpp/mersennetwister.cpp
)

find_library(FOLLY_LIBRARIES
  NAMES jemalloc glog
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
    $ENV{HOME}/.linuxbrew/lib
)
if(FOLLY_LIBRARIES)
  string(APPEND CMAKE_CXX_FLAGS " -D_USE_FOLLY")
  add_library(fbvector STATIC
    external-libraries/folly/folly/ScopeGuard.cpp
    external-libraries/folly/folly/lang/Assume.cpp
  )
  target_include_directories(fbvector SYSTEM PRIVATE
    external-libraries/folly
  )
  target_link_libraries(fbvector PRIVATE
    -lglog
    -ljemalloc
  )
  set(METALIBS ${METALIBS} fbvector)
else()
  message(WARNING "Folly (Facebook Open Source Library) not found."
  " Install folly to reduce memory taken up by the Row-Compressed annotator"
  " and enable other optimizations.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/external-libraries/binrel_wt/cmake/modules"
)
set(BRWT_INCLUDE_TESTS FALSE)
add_subdirectory(
  external-libraries/binrel_wt
)

message(STATUS "Compile flags used:"
  "${CMAKE_CXX_FLAGS}"
  " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
)

file(GLOB metagraph_core_files
  "src/*.cpp"
  "src/*/*.cpp"
  "src/*/*/*.cpp"
  "src/*/*/*/*.cpp"
)

# for compiling faster
list(REVERSE metagraph_core_files)

list(FILTER metagraph_core_files EXCLUDE REGEX ".*\\._.*")
list(FILTER metagraph_core_files EXCLUDE REGEX ".*main.cpp")
list(FILTER metagraph_core_files EXCLUDE REGEX ".*/cli/.*")

add_library(metagraph-core STATIC ${metagraph_core_files})

file(GLOB metagraph_cli_files
  "src/cli/*.cpp"
  "src/cli/*/*.cpp"
  "src/cli/*/*/*.cpp"
)
add_library(metagraph-cli STATIC ${metagraph_cli_files})

add_executable(metagraph "src/main.cpp")
set_target_properties(metagraph PROPERTIES OUTPUT_NAME "metagraph_${CMAKE_DBG_ALPHABET}")

set(METALIBS ${METALIBS}
  -lKMC
  -lhts -lz -lbz2
  -lsdsl
  brwt
  jsoncpp_lib
  progress_bar
  fmt::fmt
  mersenne_twister
  spdlog::spdlog
)
if(APPLE)
  set(METALIBS ${METALIBS} ${OpenMP_CXX_LIBRARIES})
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(METALIBS ${METALIBS} -latomic)
endif()

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set(OPENSSL_USE_STATIC_LIBS TRUE)
  find_package(OpenSSL REQUIRED)
  set(METALIBS ${METALIBS} -lbz2 -llzma -lcurl -lssl -lcrypto)
endif()

# check for std::filesystem::temp_directory_path
include(CheckCXXSourceRuns)
set(CMAKE_REQUIRED_FLAGS " -std=c++17")
check_cxx_source_runs("
    #include <iostream>
    #include <filesystem>
    int main() {
        std::cout << std::filesystem::temp_directory_path();
        return 0;
    }
" CPPNOFS)
if(NOT CPPNOFS)
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
endif()

target_link_libraries(metagraph-core ${METALIBS})
target_link_libraries(metagraph-cli PRIVATE metagraph-core)
target_link_libraries(metagraph metagraph-core metagraph-cli)

add_custom_target(link_target ALL
  COMMAND ${CMAKE_COMMAND} -E
  create_symlink
    "metagraph_${CMAKE_DBG_ALPHABET}"
    "metagraph"
)


#-------------------
# Experiments
#-------------------
file(GLOB experiment_files
  "experiments/*.cpp"
)
list(FILTER experiment_files EXCLUDE REGEX ".*\\._.*")

include_directories(
  experiments
  external-libraries/tclap/include
)

link_directories(
  external-libraries/sdsl-lite/lib
)
add_executable(run_experiments ${experiment_files})

target_link_libraries(run_experiments metagraph-core metagraph-cli)

#-------------------
# install library
#-------------------

file(GLOB header_files
  "src/*.hpp"
  "src/*/*.hpp"
  "src/*/*/*.hpp"
)

install(TARGETS metagraph-core DESTINATION lib)
install(FILES ${header_files} DESTINATION include/metagraph)
install(TARGETS metagraph DESTINATION bin)

#-------------------
# Unit Tests
#-------------------
# There is an old gtest library shipped with sdsl.
# Remove it as we are going to use a more recent one.
execute_process(COMMAND rm -rf ${PROJECT_SOURCE_DIR}/external-libraries/sdsl-lite/include/gtest)

add_subdirectory(${PROJECT_SOURCE_DIR}/external-libraries/googletest EXCLUDE_FROM_ALL)
target_compile_options(gtest_main PRIVATE -w)
target_compile_options(gtest PRIVATE -w)

set(DEATH_TEST_FLAG "")
if(NO_DEATH_TESTS)
  message(STATUS "Disabling death tests in gtest")
  set(DEATH_TEST_FLAG -D_NO_DEATH_TEST)
endif()

file(GLOB unit_tests_files
  "tests/*.cpp"
  "tests/*/*.cpp"
  "tests/*/*/*.cpp"
  "tests/*/*/*/*.cpp"
)
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files})
target_include_directories(unit_tests PRIVATE "${PROJECT_SOURCE_DIR}")
target_compile_definitions(unit_tests PRIVATE TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")
target_link_libraries(unit_tests gtest_main gtest gmock metagraph-core)

target_compile_options(unit_tests PRIVATE -Wno-undefined-var-template ${DEATH_TEST_FLAG})

#-------------------
# Benchmarks
#-------------------
add_subdirectory(${PROJECT_SOURCE_DIR}/external-libraries/benchmark EXCLUDE_FROM_ALL)

file(GLOB benchmark_files
  "benchmarks/*.cpp"
  "benchmarks/*/*.cpp"
)
list(FILTER benchmark_files EXCLUDE REGEX ".*\\._.*")

add_executable(benchmarks ${benchmark_files})
target_compile_definitions(benchmarks PRIVATE TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")
target_link_libraries(benchmarks benchmark_main benchmark metagraph-core metagraph-cli)

target_compile_options(benchmarks PRIVATE -Wno-undefined-var-template)

#-------------------
# Tests
#-------------------
find_package(Python3 3.4 REQUIRED COMPONENTS Interpreter)

add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/test_venv/DONE
  COMMAND ${PROJECT_SOURCE_DIR}/integration_tests/install_test_environment.sh ${PROJECT_SOURCE_DIR}/api/python ${PROJECT_BINARY_DIR}/test_venv
)
execute_process(COMMAND ${PROJECT_SOURCE_DIR}/integration_tests/install_test_environment.sh ${PROJECT_SOURCE_DIR}/api/python ${PROJECT_BINARY_DIR}/test_venv)

# creating a make target, s.t. the environment can be set up explicitly
add_custom_target(install_test_environment
  DEPENDS ${PROJECT_BINARY_DIR}/test_venv/DONE
)

# --- integration tests ---
add_custom_target(check
  COMMAND ${PROJECT_BINARY_DIR}/test_venv/bin/python ${PROJECT_SOURCE_DIR}/integration_tests/main.py
  DEPENDS metagraph
  DEPENDS install_test_environment
)
# --- integration tests + unit tests ---
add_custom_target(check_all
  COMMAND make check
  COMMAND unit_tests
  DEPENDS check unit_tests
)
add_test(NAME integration_tests COMMAND make check)
add_test(NAME unit_tests COMMAND unit_tests)

# --- integration tests runner ---
get_target_property(PythonPath Python3::Interpreter IMPORTED_LOCATION)
file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests" "#!/bin/bash\n")
file(APPEND "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests"
  "\n${PROJECT_BINARY_DIR}/test_venv/bin/python ${PROJECT_SOURCE_DIR}/integration_tests/main.py \"$@\""
)
file(INSTALL "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests"
  DESTINATION "${PROJECT_BINARY_DIR}/"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                   GROUP_READ GROUP_EXECUTE
                   WORLD_READ WORLD_EXECUTE
)

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
  SWIG_LINK_LIBRARIES(metagraph metagraph-core ${METALIBS} ${PYTHON_LIBRARIES})
endif()
