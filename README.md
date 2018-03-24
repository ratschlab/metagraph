# Metagenome Graph Project


### Prerequisites
- cmake 3.6.1
- C++14
- HTSlib

Can be installed with `brew` or `linuxbrew`.

### Install
1. `git clone --recursive https://github.com/ratschlab/projects2014-metagenome.git`
2. install **libmaus2** and **sdsl-lite** in `metagraph/external-libraries/` following the corresponding istructions  
or simply run the following script from the `metagraph/` dir
```bash
cd external-libraries/sdsl-lite
./install.sh $(pwd)

cd ../libmaus2
libtoolize
aclocal
autoreconf -i -f
./configure --prefix=$(pwd)
make -j $(($(getconf _NPROCESSORS_ONLN) - 1))
make install
cd ../../../
```
use `glibtoolize` instead of `libtoolize` on MacOS

3. go to the **build** directory `mkdir -p metagraph/build && cd metagraph/build`
4. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))`
5. run unit tests `./unit_tests`

### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (`Debug` by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)
- `-DPYTHON_INTERFACE=ON` -- compile python interface (requires shared libraries, OFF by default)

### Compile in the amino acid mode
- `-DCMAKE_DBG_ALPHABET=Protein` -- use the alphabet of amino acids (`DNA` by default)
