# Metagenome Graph Project


Prerequisites
- C++ 14
- boost
- cmake 3.6.1

Install
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions
3. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
4. compile by `cmake .. && make -j30 && ./unit_tests`
