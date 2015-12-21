SEQAN_DIR=/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include
LIBMAUS_DIR=/cbio/grlab/home/akahles/software/libmaus2-lib

CXX=g++ -std=c++11 -fopenmp -O3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(SEQAN_DIR) -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib -L$(LIBMAUS_DIR)/lib/libmaus2/2.0.33
LDFLAGS=-lrt -lmaus2 

default: all

all: metagraph

debug: metagraph_debug

metagraph: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp $(LDFLAGS) -D DBGDEBUG

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
