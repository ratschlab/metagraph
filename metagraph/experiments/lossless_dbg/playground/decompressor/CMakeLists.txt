cmake_minimum_required(VERSION 3.0)
project(decompressor)
add_executable(decompressor main.cpp)
target_link_libraries(
        decompressor metagraph-core ${METALIBS}
)
