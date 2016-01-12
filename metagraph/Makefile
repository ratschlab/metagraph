SEQAN_DIR=external-libraries/seqan-include
LIBMAUS_DIR=external-libraries/libmaus2-lib

CXX=g++ -std=c++11 -fopenmp -O3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(SEQAN_DIR) -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib
LDFLAGS=-lrt -lmaus2

default: all

all: metagraph

debug: metagraph_debug

metagraph: dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp metagraph.cpp 
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp $(LDFLAGS) -D DBGDEBUG

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
