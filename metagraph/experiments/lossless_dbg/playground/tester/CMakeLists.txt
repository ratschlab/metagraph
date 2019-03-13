cmake_minimum_required(VERSION 3.0)
project(sampler)
add_executable(tester main.cpp)
target_link_libraries(
        tester metagraph ${METALIBS}
)
