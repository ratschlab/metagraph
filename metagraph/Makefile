LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB_DIR=external-libraries/rocksdb
SDSL_DIR=external-libraries/local

CXX=g++ -std=c++11 -fopenmp -D_THREAD_SAFE -std=c++0x -pthread -O3 -g -msse2  #-march=core-avx-i -mtune=core-avx-i -mavx2 #-gdwarf-3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib -I$(SDSL_DIR)/include -L$(SDSL_DIR)/lib #-I$(ROCKSDB_DIR)/include -L$(ROCKSDB_DIR)
LDFLAGS=-lrt -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256 -lsdsl -lhts#-lrocksdb -lsnappy -lbz2

default: all

all: metagraph

debug: metagraph_debug

metagraph: dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp dbg_succinct_boost.hpp config.hpp metagraph.cpp config.cpp serialization.hpp datatypes.hpp helpers.hpp annotation.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp dbg_succinct_boost.hpp config.hpp serialization.hpp datatypes.hpp helpers.hpp annotation.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS) -D DBGDEBUG

bitest: bitest.cpp
	g++ -g -O3 -D__STDC_LIMIT_MACROS bitest.cpp -lhts -lz

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
