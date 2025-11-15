cmake_minimum_required(VERSION 3.16...3.30)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
cmake_policy(SET CMP0015 NEW)

# Read version from JSON
# Use CMAKE_CURRENT_SOURCE_DIR to work correctly when used as a submodule
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/../package.json MY_JSON_STRING)
string(JSON VERSION GET ${MY_JSON_STRING} version)
message(STATUS "Current MetaGraph version: ${VERSION}")

project(MetaGraph
  VERSION ${VERSION}
  DESCRIPTION "A tool for scalable construction of annotated genome graphs and sequence-to-graph alignment"
)

set(METAGRAPH_PUBLIC_DEFS)
list(APPEND METAGRAPH_PUBLIC_DEFS VERSION="${VERSION}")
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

# Option to build the metagraph executable and tests
# Default is ON when MetaGraph is the top-level project, OFF otherwise
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  option(METAGRAPH_BUILD_EXECUTABLE "Build the metagraph CLI executable" ON)
else()
  option(METAGRAPH_BUILD_EXECUTABLE "Build the metagraph CLI executable" OFF)
endif()

if(BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(OPENSSL_USE_STATIC_LIBS TRUE)
endif()

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

if(NOT CMAKE_DBG_ALPHABET)
  set(CMAKE_DBG_ALPHABET "DNA")
endif()

if(CMAKE_DBG_ALPHABET STREQUAL "DNA")
  set(METAGRAPH_ALPHABET_DEF _DNA_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA5")
  set(METAGRAPH_ALPHABET_DEF _DNA5_GRAPH)
elseif(CMAKE_DBG_ALPHABET STREQUAL "DNA_CASE_SENSITIVE")
  set(METAGRAPH_ALPHABET_DEF _DNA_CASE_SENSITIVE_GRAPH)
else()
  set(METAGRAPH_ALPHABET_DEF _PROTEIN_GRAPH)
endif()
message(STATUS "Compile for ${BoldCyan}${CMAKE_DBG_ALPHABET}${ColourReset} alphabet")
list(APPEND METAGRAPH_PUBLIC_DEFS ${METAGRAPH_ALPHABET_DEF})

set(JSONCPP_WITH_TESTS FALSE)
add_subdirectory(external-libraries/jsoncpp)

add_compile_options(
  -Wall -Wextra -Werror
  -D_THREAD_SAFE -pthread
  # -DDBGDEBUG -Wconversion -Wsign-conversion -Werror=shadow=compatible-local
)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(
    -Wno-exit-time-destructors
    -Wno-deprecated-declarations
    -Wno-vla-extension
    -Wno-deprecated-copy
  )

  if (NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    add_compile_options(
     -Wno-suggest-destructor-override
     -Wno-deprecated-declarations
    )
  endif()
endif()

if(APPLE)
  # add path to system libraries but avoid cross-compilation errors
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-Wno-poison-system-directories SUPPORT_NO_POISON_SYS_DIR)
  if(SUPPORT_NO_POISON_SYS_DIR)
    add_compile_options(-Wno-poison-system-directories)
  endif()
endif()

findOpenMP()

if(APPLE)
  if(CMAKE_OSX_ARCHITECTURES MATCHES "")
    set(CMAKE_OSX_ARCHITECTURES ${CMAKE_SYSTEM_PROCESSOR})
  elseif(${CMAKE_OSX_ARCHITECTURES} MATCHES "arm64")
    set(CMAKE_SYSTEM_PROCESSOR arm64)
    set(CMAKE_HOST_SYSTEM_PROCESSOR arm64)
  endif()
endif()

# Check if 128 bit integers are supported
try_compile(MODE_TI
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/external-libraries/sdsl-lite/CMakeModules/check_mode_ti.cpp"
)
if(MODE_TI)
  list(APPEND METAGRAPH_PUBLIC_DEFS MODE_TI)
  if(NOT ${CMAKE_SYSTEM_PROCESSOR} MATCHES "^(arm|aarch64|arm64)")
      add_compile_options(-mcx16)
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND NOT ${CMAKE_OSX_ARCHITECTURES} MATCHES "arm64")
  execute_process(COMMAND sysctl -n machdep.cpu.features
                  OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND sysctl -n machdep.cpu.leaf7_features
                  OUTPUT_VARIABLE _cpuinfo2 OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(_cpuinfo "${_cpuinfo} ${_cpuinfo2}")
endif()

if(WITH_MSSE42 OR (NOT DEFINED WITH_MSSE42 AND _cpuinfo MATCHES "(sse4_2)|(sse4a)|(SSE4.2)"))
add_compile_options(-msse4.2)
else()
message(WARNING "Compiling without msse4.2 instructions!")
endif()

if(WITH_AVX OR (NOT DEFINED WITH_AVX AND _cpuinfo MATCHES "(avx2)|(AVX2)"))
add_compile_options(-mavx -mavx2 -mfma -mbmi -mbmi2)
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
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/KMC
)
if(result)
  message(FATAL_ERROR "CMake step for KMC failed: ${result}")
endif()
include(ProcessorCount)
ProcessorCount(N)
if(N EQUAL 0)
  set(N 1)
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release -- -j${N}
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/KMC
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

# Build sdsl-lite out-of-tree and install locally to avoid target clashes
include(ExternalProject)
set(SDSL_INSTALL ${PROJECT_BINARY_DIR}/external-libraries/sdsl)
ExternalProject_Add(SDSL
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/sdsl-lite
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL}
  INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
          COMMAND ${CMAKE_COMMAND} -E rm -rf ${SDSL_INSTALL}/include/gtest
)


add_library(ips4o INTERFACE)
target_include_directories(ips4o INTERFACE external-libraries/ips4o)

add_library(caches INTERFACE)
target_include_directories(caches INTERFACE external-libraries/caches/include)

add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE external-libraries/eigen)
target_compile_options(eigen INTERFACE -Wno-unused-but-set-variable -Wno-unused-parameter)

# Allow MetaGraph's custom Find*.cmake modules to be picked up even when the
# project is included as a subdirectory.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(zstd)
if(ZSTD_FOUND)
  list(APPEND METAGRAPH_PUBLIC_DEFS _SUPPORT_ZSTD)
else()
  message(WARNING "Disabling zstd support.")
endif()

if(APPLE)
  if(NOT DEFINED ICU_ROOT)
    if(EXISTS /opt/homebrew/opt/icu4c)
      set(ICU_ROOT "/opt/homebrew/opt/icu4c")
    elseif(EXISTS /usr/local/opt/icu4c)
      set(ICU_ROOT "/usr/local/opt/icu4c")
    else()
      find_program(CONDA_EXECUTABLE conda)
      if(CONDA_EXECUTABLE)
        execute_process(COMMAND ${CONDA_EXECUTABLE} info --base
                        OUTPUT_VARIABLE ICU_ROOT
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
      endif()
    endif()
  endif()
  find_package(ICU COMPONENTS data REQUIRED)
  get_filename_component(ICU_DATA_LIBRARY_PATH ${ICU_DATA_LIBRARY} DIRECTORY)
endif()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS iostreams)

# makes ASIO not depend on boost (used by Simple-Web-Server via CMake variable)
set(USE_STANDALONE_ASIO ON)
set(ASIO_PATH external-libraries/asio/asio)

# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()
add_subdirectory(external-libraries/spdlog)
add_subdirectory(external-libraries/sshash SYSTEM)
add_subdirectory(external-libraries/DYNAMIC)
add_subdirectory(external-libraries/zlib)
target_compile_options(zlib
  PRIVATE
    -Wno-unused-variable
    -Wno-sign-compare
    -Wno-implicit-fallthrough
    -Wno-deprecated-non-prototype
    -Wno-unknown-warning-option
    -Wno-unused-parameter
    -Wno-maybe-uninitialized
)

add_library(folly STATIC
  external-libraries/folly/folly/lang/SafeAssert.cpp
  external-libraries/folly/folly/lang/ToAscii.cpp
)
target_include_directories(folly PRIVATE
  external-libraries/folly
  ${PROJECT_BINARY_DIR}/external-libraries/folly/include
)

add_library(progress_bar STATIC
  external-libraries/cpp_progress_bar/progress_bar.cpp
)

add_library(mersenne_twister STATIC
  external-libraries/rollinghashcpp/mersennetwister.cpp
)

add_library(sdust STATIC
  external-libraries/sdust/sdust.c
)

target_compile_options(sdust PRIVATE -Wall -Wc++-compat -Wno-error)

# (macOS) make arch flags string like "-arch arm64" or "-arch arm64 -arch x86_64"
function(make_osx_arch_flags outvar)
  if(APPLE)
    message(STATUS "Compiling for architecture: ${CMAKE_OSX_ARCHITECTURES}")
    set(f "")
    foreach(a ${CMAKE_OSX_ARCHITECTURES})
      list(APPEND f -arch ${a})
    endforeach()
    string(REPLACE ";" " " f "${f}")
    set(${outvar} "${f}" PARENT_SCOPE)
  else()
    set(${outvar} "" PARENT_SCOPE)
  endif()
endfunction()
make_osx_arch_flags(OSX_ARCH_FLAGS)

set(HTSLIB_INSTALL ${PROJECT_BINARY_DIR}/external-libraries/htslib)
ExternalProject_Add(HTSLIB
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/htslib
  PREFIX "${PROJECT_BINARY_DIR}/HTSLIB-build"
  CONFIGURE_COMMAND
    ${CMAKE_COMMAND} -E env
      "CPPFLAGS=${OSX_ARCH_FLAGS}\ -I${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/zlib"
      "CFLAGS=${OSX_ARCH_FLAGS}\ -I${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/zlib"
      "LDFLAGS=${OSX_ARCH_FLAGS}\ -L${PROJECT_BINARY_DIR}/external-libraries/zlib"
    /bin/sh -c "autoreconf -i ${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/htslib && exec ${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/htslib/configure --disable-libcurl --disable-bz2 --disable-lzma"
  BUILD_COMMAND $(MAKE) -j${N} lib-static
  #INSTALL_COMMAND $(MAKE) install prefix=${HTSLIB_INSTALL}  # this installs both static and dynamic
  # Install only static lib to avoid the dependency
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E make_directory       ${HTSLIB_INSTALL}/lib
    && ${CMAKE_COMMAND} -E make_directory    ${HTSLIB_INSTALL}/include/htslib
    && ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/HTSLIB-build/src/HTSLIB-build/libhts.a ${HTSLIB_INSTALL}/lib/
    && ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/external-libraries/htslib/htslib ${HTSLIB_INSTALL}/include/htslib
  DEPENDS zlib
)

set(OTHER_COMPILER_FLAGS "")
get_directory_property(DIRECTORY_DEFINITIONS COMPILE_DEFINITIONS)
foreach(DEF ${DIRECTORY_DEFINITIONS})
  if(DEF)
    string(REPLACE "\"" "\\\"" DEF "${DEF}")
    list(APPEND OTHER_COMPILER_FLAGS " -D${DEF}")
  endif()
endforeach()
get_directory_property(DIRECTORY_OPTIONS COMPILE_OPTIONS)
foreach(OPT ${DIRECTORY_OPTIONS})
  if(OPT)
    string(REPLACE "\"" "\\\"" OPT "${OPT}")
    list(APPEND OTHER_COMPILER_FLAGS " ${OPT}")
  endif()
endforeach()
message(STATUS "Compile flags used:"
  "${CMAKE_CXX_FLAGS}"
  " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
  " ${OTHER_COMPILER_FLAGS}"
)

# Create metagraph-core target early so we can use target_include_directories
file(GLOB_RECURSE metagraph_core_files "src/*.cpp")
list(REVERSE metagraph_core_files)  # for compiling faster
list(FILTER metagraph_core_files EXCLUDE REGEX ".*\\._.*")
list(FILTER metagraph_core_files EXCLUDE REGEX ".*main.cpp")
list(FILTER metagraph_core_files EXCLUDE REGEX ".*/cli/.*")
add_library(metagraph-core STATIC ${metagraph_core_files})

# Set up include directories for metagraph-core
target_include_directories(metagraph-core PUBLIC
  ${PROJECT_SOURCE_DIR}/src
  external-libraries/KMC/kmc_api
  ${SDSL_INSTALL}/include
  external-libraries/hopscotch-map/include
  external-libraries/ordered-map/include
  external-libraries/cpp_progress_bar
  external-libraries/Priority-Deque
  external-libraries/rollinghashcpp
  external-libraries/spdlog/include
  ${PROJECT_BINARY_DIR}/external-libraries/htslib/include
  external-libraries/zlib
  external-libraries/sdust
  external-libraries/jsoncpp/include
  external-libraries/DYNAMIC/include
)

# Propagate OpenMP flags to consumers (compile and link)
# GCC needs -fopenmp during both compilation and linking to auto-link -lgomp
if(OpenMP_CXX_FLAGS)
  # Convert space-separated string to list for proper handling by target_*_options
  separate_arguments(OpenMP_CXX_FLAGS_LIST UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
  target_compile_options(metagraph-core PUBLIC ${OpenMP_CXX_FLAGS_LIST})
  target_link_options(metagraph-core PUBLIC ${OpenMP_CXX_FLAGS_LIST})
endif()

target_link_directories(metagraph-core INTERFACE
  ${PROJECT_BINARY_DIR}/KMC
  ${HTSLIB_INSTALL}/lib
  ${SDSL_INSTALL}/lib
)

if(APPLE)
  target_link_directories(metagraph-core INTERFACE /usr/local/lib)
  if(EXISTS /opt/homebrew/lib)
    target_link_directories(metagraph-core INTERFACE /opt/homebrew/lib)
  endif()
  target_include_directories(metagraph-core SYSTEM PUBLIC
    /usr/local/include
    /opt/homebrew/include
  )
endif()

if(ICU_DATA_LIBRARY_PATH)
  target_link_directories(metagraph-core INTERFACE ${ICU_DATA_LIBRARY_PATH})
endif()

# TODO: move ASIO out of here once it replaces sprintf with snprintf
target_include_directories(metagraph-core SYSTEM PUBLIC
  external-libraries/asio/asio/include
  external-libraries/simde-no-tests
)

# Configure Folly if available
findFolly()
if(FOLLY_FOUND)
  list(APPEND METAGRAPH_PUBLIC_DEFS _USE_FOLLY)
  if(Jemalloc_FOUND)
    list(APPEND METAGRAPH_PUBLIC_DEFS USE_JEMALLOC)
  else()
    message(WARNING "Jemalloc not found, some optimizations will be disabled.")
  endif()
  set(METALIBS ${METALIBS} fbvector)
  target_include_directories(metagraph-core SYSTEM PUBLIC
    external-libraries/folly
    ${PROJECT_BINARY_DIR}/external-libraries/folly/include
  )
endif()

# Propagate public compile definitions to dependents during build
target_compile_definitions(metagraph-core PUBLIC ${METAGRAPH_PUBLIC_DEFS})

# Generate config header with all compile definitions for external consumers
# External consumers will get definitions from the header instead of CMake
# Convert list items to #define directives - process both the variable and target property
get_target_property(TARGET_DEFS metagraph-core COMPILE_DEFINITIONS)
if(TARGET_DEFS)
  list(APPEND METAGRAPH_PUBLIC_DEFS ${TARGET_DEFS})
  list(REMOVE_DUPLICATES METAGRAPH_PUBLIC_DEFS)
endif()

foreach(def ${METAGRAPH_PUBLIC_DEFS})
  if(def MATCHES "^([^=]+)=(.*)$")
    # Definition with value: FOO=bar -> #define FOO bar
    set(${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
  else()
    # Boolean definition: FOO -> #define FOO
    set(${def} 1)
  endif()
endforeach()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/common/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/gen/common/config.hpp
    @ONLY
)

# Add the generated header to include path
target_include_directories(metagraph-core PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/gen)

file(GLOB_RECURSE metagraph_cli_files "src/cli/*.cpp")
add_library(metagraph-cli STATIC ${metagraph_cli_files})
add_dependencies(metagraph-core HTSLIB SDSL)


set(METALIBS ${METALIBS}
  -lKMC
  -lhts
  zlib
  deflate
  Boost::iostreams
  -lsdsl
  jsoncpp_lib
  progress_bar
  mersenne_twister
  sdust
  spdlog::spdlog
  sshash_static
  XXSDS_DYNAMIC
  ips4o
  caches
  eigen
  folly
)
if(APPLE)
  set(METALIBS ${METALIBS} ${OpenMP_CXX_LIBRARIES})
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(METALIBS ${METALIBS} -latomic)
endif()

if(ZSTD_FOUND)
  set(METALIBS ${METALIBS} zstd::zstd)
endif()

if(BUILD_STATIC AND NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libstdc++ -static-libgcc")
  find_package(OpenSSL REQUIRED)
  set(METALIBS ${METALIBS} -lssl -lcrypto -lgomp -lrt -ldl -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
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

target_link_libraries(metagraph-core PUBLIC ${METALIBS})
add_subdirectory(external-libraries/Simple-Web-Server)
target_link_libraries(metagraph-cli PUBLIC metagraph-core PRIVATE simple-web-server)

if(METAGRAPH_BUILD_EXECUTABLE)
  add_executable(metagraph "src/main.cpp")
  set_target_properties(metagraph PROPERTIES OUTPUT_NAME "metagraph_${CMAKE_DBG_ALPHABET}")
  target_link_libraries(metagraph metagraph-core metagraph-cli)

  add_custom_target(link_target ALL
    COMMAND ${CMAKE_COMMAND} -E
    create_symlink
      "metagraph_${CMAKE_DBG_ALPHABET}"
      "metagraph"
  )

  add_dependencies(metagraph link_target)
endif()


#-------------------
# install library
#-------------------

get_target_property(INCLUDE_DIRS metagraph-core INCLUDE_DIRECTORIES)

# Run the header installation at install time (after ExternalProjects are built)
install(CODE "
  set(INCLUDE_DIRS \"${INCLUDE_DIRS}\")
  set(PROJECT_SOURCE_DIR \"${PROJECT_SOURCE_DIR}\")
  
  foreach(dir \${INCLUDE_DIRS})
    if(\"\${dir}\" STREQUAL \"\${PROJECT_SOURCE_DIR}\")
      continue()
    endif()
    
    if(NOT EXISTS \"\${dir}\")
      continue()
    endif()
    
    file(GLOB_RECURSE subdirs LIST_DIRECTORIES true \"\${dir}/*\")
    list(APPEND subdirs \${dir})
    
    foreach(subdir \${subdirs})
      if(NOT IS_DIRECTORY \"\${subdir}\")
        continue()
      endif()
      
      file(GLOB headers \"\${subdir}/*.h\" \"\${subdir}/*.hpp\")
      if(headers)
        file(RELATIVE_PATH rel_path \"\${dir}\" \"\${subdir}\")
        file(INSTALL \${headers}
             DESTINATION \"\${CMAKE_INSTALL_PREFIX}/include/\${rel_path}\")
      endif()
    endforeach()
  endforeach()
")

# Collect static libraries from metagraph-core's dependencies
# Use $<TARGET_FILE:...> to get actual library paths
set(BUNDLE_LIBS_LIST "$<TARGET_FILE:metagraph-core>")

get_target_property(LINK_LIBS metagraph-core LINK_LIBRARIES)
if(LINK_LIBS)
  foreach(LIB ${LINK_LIBS})
    # Skip if not a CMake target
    if(NOT TARGET ${LIB})
      continue()
    endif()
    
    # Check if it's a static library target
    get_target_property(LIB_TYPE ${LIB} TYPE)
    if(LIB_TYPE STREQUAL "STATIC_LIBRARY")
      # Skip imported libraries (e.g., system Boost)
      get_target_property(IS_IMPORTED ${LIB} IMPORTED)
      if(IS_IMPORTED)
        continue()
      endif()
      
      # Add our internal static library
      string(APPEND BUNDLE_LIBS_LIST ";$<TARGET_FILE:${LIB}>")
    endif()
  endforeach()
endif()

# Add external libraries built by ExternalProject (not CMake targets)
string(APPEND BUNDLE_LIBS_LIST ";${SDSL_INSTALL}/lib/libsdsl.a")
string(APPEND BUNDLE_LIBS_LIST ";${HTSLIB_INSTALL}/lib/libhts.a")
string(APPEND BUNDLE_LIBS_LIST ";${PROJECT_BINARY_DIR}/KMC/libKMC.a")

# Bundle all internal dependencies into libmetagraph-core.a during install
# This creates a combined static library containing metagraph-core and all bundled dependencies
install(CODE "
  set(PROJECT_BINARY_DIR \"${PROJECT_BINARY_DIR}\")
  set(PROJECT_SOURCE_DIR \"${PROJECT_SOURCE_DIR}\")
  set(SDSL_INSTALL \"${SDSL_INSTALL}\")
  set(HTSLIB_INSTALL \"${HTSLIB_INSTALL}\")
  set(CMAKE_AR \"${CMAKE_AR}\")
  set(BUNDLE_LIBS_LIST \"${BUNDLE_LIBS_LIST}\")
  
  message(STATUS \"Bundling dependencies into libmetagraph-core-bundle.a...\")
  
  # Create temporary directory for bundling
  file(REMOVE_RECURSE \${PROJECT_BINARY_DIR}/metagraph-core-bundle)
  file(MAKE_DIRECTORY \${PROJECT_BINARY_DIR}/metagraph-core-bundle)
  
  # Extract each library into its own subdirectory to avoid filename collisions
  foreach(LIB \${BUNDLE_LIBS_LIST})
    if(NOT EXISTS \"\${LIB}\")
      message(FATAL_ERROR \"Required library not found for bundling: \${LIB}\")
    endif()
    message(STATUS \"  Extracting: \${LIB}\")
    
    # Get library name without path and extension for subdirectory
    get_filename_component(LIB_NAME \"\${LIB}\" NAME_WE)
    string(REGEX REPLACE \"^lib\" \"\" LIB_NAME \"\${LIB_NAME}\")
    
    # Create subdirectory for this library
    file(MAKE_DIRECTORY \${PROJECT_BINARY_DIR}/metagraph-core-bundle/\${LIB_NAME})
    
    # Extract into subdirectory
    execute_process(COMMAND \${CMAKE_AR} x \${LIB}
                    WORKING_DIRECTORY \${PROJECT_BINARY_DIR}/metagraph-core-bundle/\${LIB_NAME}
                    RESULT_VARIABLE AR_RESULT)
    if(AR_RESULT)
      message(FATAL_ERROR \"Failed to extract \${LIB} (exit code: \${AR_RESULT})\")
    endif()
  endforeach()
  
  # Create the combined archive with -bundle suffix to avoid overwriting the original
  execute_process(COMMAND /bin/bash \${PROJECT_SOURCE_DIR}/bundle_static_libs.sh 
                          \${PROJECT_BINARY_DIR}/metagraph-core-bundle 
                          \${PROJECT_BINARY_DIR}/libmetagraph-core-bundle.a 
                          \${CMAKE_AR}
                  RESULT_VARIABLE BUNDLE_RESULT)
  if(BUNDLE_RESULT)
    message(FATAL_ERROR \"Failed to bundle static libraries (exit code: \${BUNDLE_RESULT})\")
  endif()
  
  if(NOT EXISTS \"\${PROJECT_BINARY_DIR}/libmetagraph-core-bundle.a\")
    message(FATAL_ERROR \"Bundled libmetagraph-core-bundle.a was not created\")
  endif()
  
  execute_process(COMMAND du -h \${PROJECT_BINARY_DIR}/libmetagraph-core-bundle.a
                  OUTPUT_VARIABLE BUNDLE_SIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE \"^([^\\t]+).*\" \"\\\\1\" BUNDLE_SIZE \"\${BUNDLE_SIZE}\")
  message(STATUS \"Created libmetagraph-core-bundle.a (\${BUNDLE_SIZE})\")
")

# Install the bundled library as libmetagraph-core.a (without -bundle suffix)
install(FILES ${PROJECT_BINARY_DIR}/libmetagraph-core-bundle.a DESTINATION lib RENAME libmetagraph-core.a)
install(TARGETS metagraph-cli DESTINATION lib)

if(METAGRAPH_BUILD_EXECUTABLE)
  install(TARGETS metagraph DESTINATION bin)

  #-------------------
  # Unit Tests
  #-------------------

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

  add_executable(unit_tests ${unit_tests_files})
  target_link_directories(unit_tests PRIVATE ${PROJECT_BINARY_DIR})
  target_include_directories(unit_tests PRIVATE "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(unit_tests PRIVATE TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")
  target_link_libraries(unit_tests gtest_main gtest gmock metagraph-core metagraph-cli)

  target_compile_options(unit_tests PRIVATE -Wno-uninitialized ${DEATH_TEST_FLAG})

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
  target_link_libraries(benchmarks benchmark_main benchmark metagraph-core metagraph-cli -lsdsl)

  target_compile_options(benchmarks PRIVATE)

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
  add_custom_target(integration_tests
    COMMAND ${PROJECT_BINARY_DIR}/test_venv/bin/python ${PROJECT_SOURCE_DIR}/integration_tests/main.py
    DEPENDS metagraph
    DEPENDS install_test_environment
  )
  # --- integration tests + unit tests ---
  add_custom_target(check_all
    COMMAND make integration_tests
    COMMAND unit_tests
    DEPENDS integration_tests unit_tests
  )
  add_test(NAME integration_tests COMMAND make integration_tests)
  add_test(NAME unit_tests COMMAND unit_tests)

  # --- integration tests runner ---
  get_target_property(PythonPath Python3::Interpreter IMPORTED_LOCATION)
  file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests" "#!/bin/bash\n")
  file(APPEND "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests"
    "\n${PROJECT_BINARY_DIR}/test_venv/bin/python ${PROJECT_SOURCE_DIR}/integration_tests/main.py \"$@\""
  )
  file(INSTALL "${PROJECT_BINARY_DIR}/CMakeFiles/integration_tests"
    DESTINATION "${PROJECT_BINARY_DIR}"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                     GROUP_READ GROUP_EXECUTE
                     WORLD_READ WORLD_EXECUTE
  )

endif() # CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME
