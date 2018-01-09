# Metagenome Graph Project


Prerequisites
- cmake 3.6.1
- C++14
- boost
- HTSlib

Install
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions
3. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
4. compile by `cmake .. && make && ./unit_tests`

Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (Debug by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)
- `-DPYTHON_INTERFACE=ON` -- compile python interface (requires shared libraries, OFF by default)

