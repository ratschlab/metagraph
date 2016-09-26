LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB_DIR=external-libraries/rocksdb

CXX=g++ -std=c++11 -fopenmp -D_THREAD_SAFE -std=c++0x -pthread -O3 #-g -gdwarf-3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib -I$(ROCKSDB_DIR)/include -L$(ROCKSDB_DIR)
LDFLAGS=-lrt -lmaus2 -lrocksdb -lsnappy -lbz2

default: all

all: metagraph

debug: metagraph_debug

metagraph: dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp metagraph.cpp config.cpp serialization.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp serialization.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS) -D DBGDEBUG

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
