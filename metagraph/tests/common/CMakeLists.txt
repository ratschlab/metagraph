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

set(EXTERNAL_LIB_DIR "../../external-libraries")
link_directories(
        ${EXTERNAL_LIB_DIR}/sdsl-lite/lib
        /usr/local/lib/
)

file(GLOB test_files "*.cpp")
add_executable(tests ${test_files})
target_include_directories(tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(tests gtest_main ${GTEST_BOTH_LIBRARIES} common)
message("Gtest libs " ${GTEST_BOTH_LIBRARIES})
