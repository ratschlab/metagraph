cmake_minimum_required(VERSION 3.0)
project(compressor)
add_executable(compressor main.cpp)
target_link_libraries(
        compressor metagraph ${METALIBS}
)
