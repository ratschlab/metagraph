cmake_minimum_required(VERSION 3.6.1)
project(test_common)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../../src/common test_common)
include_directories(../../src)
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

enable_testing()


file(GLOB test_files "*.cpp")
add_definitions(-D_OPENMP) # TODO(ddanciu) - this should somehow be set automatically
add_executable(test_common ${test_files})
target_include_directories(test_common SYSTEM PUBLIC "/usr/local/include")
target_link_libraries(test_common gtest_main gtest common)
