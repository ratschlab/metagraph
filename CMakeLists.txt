# Wrapper CMakeLists.txt for using metagraph as a submodule
# This allows the metagraph/CMakeLists.txt to find package.json correctly
# when metagraph is included via add_subdirectory() from another project

cmake_minimum_required(VERSION 3.16)

# When using metagraph as a submodule, just forward to the actual CMakeLists.txt
add_subdirectory(metagraph)
