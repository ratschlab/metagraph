CXX=g++ -std=c++11 -fopenmp -O3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
#CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-seqan-v2.0.0/include 
#CXXFLAGS+= -I/cbio/grlab/home/akahles/software/sdsl-lite-lib/include -L/cbio/grlab/home/akahles/software/sdsl-lite-lib/lib 
#CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include -I /cbio/grlab/home/akahles/software/cpp-btree-1.0.1 -I/cbio/grlab/home/akahles/software/libmaus2-lib/include -I. -L/cbio/grlab/home/akahles/software/libmaus2-lib/lib
CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include -I/cbio/grlab/home/akahles/software/libmaus2-lib/include -I. -L/cbio/grlab/home/akahles/software/libmaus2-lib/lib -L/cbio/grlab/home/akahles/software/libmaus2-lib/lib/libmaus2/2.0.33
LDFLAGS=-lrt -lmaus2 # -lsdsl -ldivsufsort -ldivsufsort64

default: all

all: metagraph

debug: metagraph_debug

metagraph: metagraph.cpp dbg_seqan.hpp dbg_succinct_libmaus.hpp dmm_tree.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_seqan.hpp dbg_succinct_libmaus.hpp dmm_tree.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp $(LDFLAGS) -D DBGDEBUG

#metagraph.o: metagraph.cpp dbg_seqan.hpp dbg_succinct_libmaus.hpp dmm_tree.hpp
#	$(CXX) $(CXXFLAGS) -c -o metagraph.o metagraph.cpp $(LDFLAGS)

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
