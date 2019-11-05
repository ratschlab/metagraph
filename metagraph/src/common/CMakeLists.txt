cmake_minimum_required(VERSION 3.6.1)
project(common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# unset OpenMP_CXX_FLAGS. If OpenMP loading fails once, it will be set to NOT_FOUND
# and happily cached, thus messing up the compilation flags forever
unset(OpenMP_CXX_FLAGS CACHE)
if ((CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
  find_package(OPENMP REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
message("OpenMp flags are " ${OpenMP_CXX_FLAGS})

if (APPLE)
  # Fix include directories on 10.14+, as they are not added by default any longer
  string(APPEND CMAKE_CXX_FLAGS " \
    -I/usr/include \
    -I/usr/local/include \
  ")
  # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
  LINK_DIRECTORIES(/usr/lib)
  LINK_DIRECTORIES(/usr/local/lib)
endif()

set(EXTERNAL_LIB_DIR "../../external-libraries")
include_directories(
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
)
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()
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
  set(METALIBS ${METALIBS}
      ${FOLLY_LIBRARIES}
      -ldouble-conversion
      -ldl
      -lboost_regex
      -lboost_filesystem
      -lboost_program_options
      -lcrypto
      ${OpenMP_CXX_LIBRARIES}
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
      $ENV{HOME}/.linuxbrew/lib
      )
  if(LIBERTY)
    set(METALIBS ${METALIBS} ${LIBERTY})
  endif()
else()
  message(WARNING "Folly (Facebook Open Source Library) not found."
      " Install folly to reduce memory taken up by the Row-Compressed annotator"
      " and enable other optimizations.")
endif()

message(STATUS "Compile flags used:"
    "${CMAKE_CXX_FLAGS}"
    " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}"
    )

file(GLOB src_files
    "*.cpp"
    "*/*.cpp"
    )

add_library(common STATIC ${src_files})
target_link_libraries(common  ${FOLLY_LIBRARIES} -lsdsl)
