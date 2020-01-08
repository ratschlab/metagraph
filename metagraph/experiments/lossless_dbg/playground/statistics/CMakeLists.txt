cmake_minimum_required(VERSION 3.0)
project(sampler)
add_executable(sampler main.cpp)
target_link_libraries(
        sampler metagraph-core ${METALIBS}
)
