cmake_minimum_required(VERSION 3.6.1)
project(metagraph_common)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    -lsdsl
    Boost::iostreams
)
