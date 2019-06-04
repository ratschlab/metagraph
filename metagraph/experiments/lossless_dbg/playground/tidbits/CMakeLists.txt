cmake_minimum_required(VERSION 3.0)
project(tidbits)
add_executable(memory_dynamic_wavelet memory_dynamic_wavelet.cpp)
target_link_libraries(memory_dynamic_wavelet metagraph ${METALIBS})
