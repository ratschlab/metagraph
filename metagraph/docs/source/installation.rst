.. _installation:

Installation
============

The core of MetaGraph is written in C++ and has been successfully tested on Linux and MacOS. In the
following, we provide detailed instructions for setting up MetaGraph.

Install with conda
------------------

There are conda packages available on bioconda for both Linux and Mac OS X::

    conda install -c bioconda -c conda-forge metagraph

The executables are called ``metagraph_DNA`` (with a ``metagraph`` symlink) and ``metagraph_Protein``.

For support of other/custom alphabets, compile from source (see `Install from source`_).


Docker container
----------------

If docker is available on your system, you can immediately get started with::

    docker run -v ${DATA_DIR_HOST}:/mnt ghcr.io/ratschlab/metagraph:master \
        build -v -k 10 -o /mnt/transcripts_1000 /mnt/transcripts_1000.fa


where ``${DATA_DIR_HOST}`` should be replaced with a directory on the host system that will be
mapped under ``/mnt`` in the container. This docker container uses the latest version of MetaGraph from
the source `GitHub repository <https://github.com/ratschlab/metagraph>`_ (branch ``master``).
See also the `image overview <https://github.com/ratschlab/metagraph/pkgs/container/metagraph>`_ for
other versions of the image.

By default, it executes the binary compiled for the DNA alphabet.
To run the binary compiled for the `Protein` alphabet, just add ``--entrypoint metagraph_Protein``::

    docker run --entrypoint metagraph_Protein \
               -v ${DATA_DIR_HOST}:/mnt ghcr.io/ratschlab/metagraph:master \
        build -v -k 10 -o /mnt/graph /mnt/protein.fa

As you can see, running MetaGraph from docker containers is very easy.
Also, the following command (or similar) may be handy to see what directory is mounted in the
container or to do other sorts of debugging::

    docker run -v ${DATA_DIR_HOST}:/mnt --entrypoint ls ghcr.io/ratschlab/metagraph:master /mnt


Install from source
-------------------

Prerequisites
^^^^^^^^^^^^^
Before compiling MetaGraph, install the following dependencies:

- cmake 3.10 or higher
- GNU GCC with C++17 (gcc-8.0.1 or higher), LLVM Clang (clang-7 or higher), or AppleClang
- bzip2

*Optional:*

- boost and jemalloc-4.0.0 or higher (to build with *folly* for efficient small vector support)
- Python 3 (for running integration tests)

.. tip:: For those without administrator/root privileges, we recommend using
         `brew <https://brew.sh/>`_ (available for MacOS and Linux).

.. tabs::

    .. group-tab:: AppleClang on MacOS

        For compiling with **AppleClang**, the prerequisites can be installed as easy as::

            brew install libomp cmake make bzip2 boost jemalloc automake autoconf


    .. group-tab:: Ubuntu / Debian

        For **Ubuntu** (20.04 LTS or higher) or **Debian** (10 or higher)::

            sudo apt-get install cmake libbz2-dev libjemalloc-dev libboost-all-dev automake autoconf


    .. group-tab:: CentOS

        For **CentOS** (8 or higher)::

            yum install cmake bzip2-devel jemalloc-devel boost-devel automake autoconf


    .. group-tab:: brew + GNU gcc

        GNU GCC and all the prerequisites can be installed with `brew <https://brew.sh/>`_ as follows::

            brew install gcc autoconf automake libtool cmake make
            [[ "$OSTYPE" == "darwin"* ]] \
                && brew remove -f boost \
                && brew install --cc=gcc-7 boost \
                && brew install gcc@9
            [[ "$OSTYPE" != "darwin"* ]] \
                && brew install gcc@9 libomp \
                && brew remove -f openssl@1.1 boost \
                && brew install --cc=gcc-9 openssl@1.1 boost

        Then, the following environment variables have to be set::

            echo "\
            # Use gcc-9 with cmake
            export CC=\"\$(which gcc-9)\"
            export CXX=\"\$(which g++-9)\"
            " >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )

    .. group-tab:: brew + LLVM Clang

        For compiling with LLVM Clang installed with `brew <https://brew.sh/>`_, the prerequisites can be installed with::

            brew install llvm libomp autoconf automake libtool cmake make boost

        Then, the following environment variables have to be set::

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
^^^^^^^^^
To compile MetaGraph, please follow these steps.

#. Clone the latest version of the code from the git repository::

    git clone --recursive https://github.com/ratschlab/metagraph.git

#. Change into the ``metagraph`` directory::

    cd metagraph

#. Make sure all submodules have been downloaded::

    git submodule update --init --recursive

#. Install *sdsl-lite* in ``metagraph/external-libraries/sdsl-lite`` with the following script::

    git submodule sync
    git submodule update --init --recursive

    pushd metagraph/external-libraries/sdsl-lite
    ./install.sh $PWD
    popd

#. Set up the ``build`` directory and change into it::

    mkdir metagraph/build
    cd metagraph/build

#. Compile::

    cmake ..
    make -j $(($(getconf _NPROCESSORS_ONLN) - 1))

#. Run unit tests (optional)::

    ./unit_tests --gtest_filter="*"

#. Run integration tests (optional)::

    ./integration_tests --test_filter="*"

Build configurations
^^^^^^^^^^^^^^^^^^^^

When configuring via ``cmake .. <arguments>`` additional arguments can be provided:

- ``-DCMAKE_BUILD_TYPE=[Debug|Release|Profile|GProfile]`` -- build modes (``Release`` by default)
- ``-DBUILD_STATIC=[ON|OFF]`` -- link statically (``OFF`` by default)
- ``-DLINK_OPT=[ON|OFF]`` -- enable link time optimization (``OFF`` by default)
- ``-DBUILD_KMC=[ON|OFF]`` -- compile the KMC executable (``ON`` by default)
- ``-DWITH_AVX=[ON|OFF]`` -- compile with *avx* instructions (``ON`` by default, if available)
- ``-DWITH_MSSE42=[ON|OFF]`` -- compile with *msse4.2* instructions (``ON`` by default, if available)
- ``-DCMAKE_DBG_ALPHABET=[Protein|DNA|DNA5|DNA_CASE_SENSITIVE]`` -- alphabet to use (``DNA`` by default)


Install API
----------------------------
See :ref:`API Install <install api>`.
