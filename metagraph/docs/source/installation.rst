.. _installation:

Installation
============

MetaGraph is written in C++ and has been successfully tested on Linux and MacOS platforms. In the
following, we provide detailed instructions for setting up the framework.

Prerequisites
-------------

Before compiling MetaGraph, you need to install the following dependencies. For users that do not
have administration/root access to their machine, we recommend the usage of `brew
<https://brew.sh/>`_ or `linuxbrew <https://linuxbrew.sh/>`_. 

- cmake 3.6.1
- GNU GCC with C++17 (gcc-9 or higher), LLVM Clang (clang-7 or higher), or AppleClang (clang-1100.0.33.8 or higher)
- HTSlib
- folly (optional)
- Python 3 (for running integration tests) 


When compiling on Mac with AppleClang
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For compiling with AppleClang, the prerequisites can be installed as easy as::

    brew install libomp cmake make htslib boost folly

When compiling on Linux with GNU GCC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For compiling with GNU GCC, the prerequisites can be installed as::

    brew install gcc autoconf automake libtool cmake make htslib
    [[ "$OSTYPE" == "darwin"* ]] \
        && brew remove -f boost double-conversion gflags glog lz4 snappy zstd folly \
        && brew install --cc=gcc-7 boost folly \
        && brew install gcc@9
    [[ "$OSTYPE" != "darwin"* ]] \
        && brew install gcc@9 libomp \
        && brew remove -f openssl@1.1 boost double-conversion gflags glog lz4 snappy zstd folly \
        && brew install --cc=gcc-5 glog zstd \
        && brew install --cc=gcc-9 openssl@1.1 boost folly

Then, the following environment variable need to be set::

    echo "\
    # Use gcc-9 with cmake
    export CC=\"\$(which gcc-9)\"
    export CXX=\"\$(which g++-9)\"
    " >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )

When compiling on Linux with LLVM Clang
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For compiling with LLVM Clang, the prerequisites can be installed as::

    brew install llvm libomp autoconf automake libtool cmake make htslib boost folly

Then, the following environment variable need to be set::

    echo "\
    # OpenMP
    export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix libomp)/lib\"
    export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix libomp)/include\"
    export CXXFLAGS=\"\$CXXFLAGS -I$(brew --prefix libomp)/include\"
    # Clang C++ flags
    export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix llvm)/lib -Wl,-rpath,$(brew --prefix llvm)/lib\"
    export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix llvm)/include\"
    export CXXFLAGS=\"\$CXXFLAGS -stdlib=libc++\"
    # Path to Clang
    export PATH=\"$(brew --prefix llvm)/bin:\$PATH\"
    # Use Clang with cmake
    export CC=\"\$(which clang)\"
    export CXX=\"\$(which clang++)\"
    " >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )

Compiling
---------
To compile MetaGraph, please follow the subsequent steps:

1. Clone the latest version of the code from the git repository::

    git clone --recursive https://github.com/ratschlab/metagraph.git

2. Make sure all submodules are downloaded:: 

    git submodule update --init --recursive

3. Install **libmaus2** and **sdsl-lite** in :code:`metagraph/external-libraries/` following the corresponding instructions or simply run the following script::

    git submodule sync
    git submodule update --init --recursive

    pushd metagraph/external-libraries/sdsl-lite
    ./install.sh $PWD
    popd

    pushd metagraph/external-libraries/libmaus2
    cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD .
    make -j $(($(getconf _NPROCESSORS_ONLN) - 1))
    make install
    popd    

4. Go to the :code:`build` directory::

    mkdir -p metagraph/build && cd metagraph/build

5. Compile::

    cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))

6. Run unit tests::

    ./unit_tests

7. Run integration tests::

    ./integration_tests

Build types
-----------

When building with :code:`cmake .. <arguments>` additional arguments can be provided as follows:

- :code:`-DCMAKE_BUILD_TYPE=[Debug|Release|Profile|GProfile]` -- build modes (:code:`Release` by default)
- :code:`-DBUILD_STATIC=[ON|OFF]` -- link statically (:code:`OFF` by default)
- :code:`-DLINK_OPT=[ON|OFF]` -- enable link time optimization (:code:`OFF` by default)
- :code:`-DBUILD_KMC=[ON|OFF]` -- compile the KMC executable (:code:`ON` by default)
- :code:`-DWITH_AVX=[ON|OFF]` -- compile with support for the avx instructions (:code:`ON` by default, if available)
- :code:`-DWITH_MSSE42=[ON|OFF]` -- compile with support for the msse4.2 instructions (:code:`ON` by default, if available)
- :code:`-DCMAKE_DBG_ALPHABET=[Protein|DNA|DNA5|DNA_CASE_SENSITIVE]` -- alphabet to use (:code:`DNA` by default)
