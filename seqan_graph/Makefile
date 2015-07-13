CXX=g++ -std=c++11 -O2
CXXFLAGS+= -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
#CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-seqan-v2.0.0/include 
CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include -I /cbio/grlab/home/akahles/software/cpp-btree-1.0.1 -I/cbio/grlab/home/akahles/software/sdsl-lite-lib/include -I. -L/cbio/grlab/home/akahles/software/sdsl-lite-lib/lib
LDFLAGS=-lrt -lsdsl -ldivsufsort -ldivsufsort64

default: all

all: metagraph

metagraph: metagraph.o dbg_seqan.hpp dbg_succinct.hpp unix_tools.hpp dmm_tree.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.o $(LDFLAGS)

metagraph.o: metagraph.cpp dbg_seqan.hpp dbg_succinct.hpp dmm_tree.hpp
	$(CXX) $(CXXFLAGS) -c -o metagraph.o metagraph.cpp $(LDFLAGS)

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
