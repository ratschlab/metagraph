cmake_minimum_required(VERSION 3.6.1)
project(metagraph_common)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED COMPONENTS iostreams)
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

set(EXTERNAL_LIB_DIR "../../external-libraries")
link_directories(
  ${EXTERNAL_LIB_DIR}/sdsl-lite/lib
)

set(FOLLY_LIBRARIES FOLLY_LIBRARIES-NOTFOUND)
find_library(FOLLY_LIBRARIES
  NAMES folly
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
    $ENV{HOME}/.linuxbrew/lib
)
if(FOLLY_LIBRARIES)
  string(APPEND CMAKE_CXX_FLAGS " -D_USE_FOLLY")
else()
  message(WARNING "Folly (Facebook Open Source Library) not found."
          " Install folly to reduce memory taken up by the Row-Compressed annotator"
          " and enable other optimizations.")
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
    set(FILESYSTEM_LIBRARY -lc++fs)
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
      set(FILESYSTEM_LIBRARY -lstdc++fs)
    else()
      message(FATAL_ERROR "std::filesystem not found")
    endif()
  endif()
endif()

file(GLOB src_files "*.cpp" "*/*.cpp")

add_library(metagraph_common STATIC ${src_files})
set(EXTERNAL_LIB_DIR "../../external-libraries")

target_include_directories(metagraph_common
  PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
    .. # this ensures that this libs include files are visible to its deps
)

target_link_libraries(metagraph_common
  PUBLIC
    ${FOLLY_LIBRARIES}
    ${FILESYSTEM_LIBRARY}
    -lsdsl
    Boost::iostreams
)
