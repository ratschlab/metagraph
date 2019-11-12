cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../../src/common test_common)

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

enable_testing()
find_package(GTest REQUIRED)

set(EXTERNAL_LIB_DIR "${CMAKE_SOURCE_DIR}/../../external-libraries")

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})
if(APPLE)
  # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
  target_link_directories(tests INTERFACE /usr/local/lib )
endif()
message("Lib dir is " ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_link_directories(tests INTERFACE ${EXTERNAL_LIB_DIR}/sdsl-lite/lib)
target_include_directories(tests
  PRIVATE
    ${GTEST_INCLUDE_DIRS}
    "../../src/"
    ${Boost_INCLUDE_DIRS}
    ${EXTERNAL_LIB_DIR}/asio/asio/include
    ${EXTERNAL_LIB_DIR}/DYNAMIC/include
    ${EXTERNAL_LIB_DIR}/hopscotch-map/include
    ${EXTERNAL_LIB_DIR}/ips4o
    ${EXTERNAL_LIB_DIR}/libmaus2/include
    ${EXTERNAL_LIB_DIR}/ordered-map/include
    ${EXTERNAL_LIB_DIR}/sdsl-lite/include
)
target_link_libraries(tests ${GTEST_BOTH_LIBRARIES} metagraph_common)
