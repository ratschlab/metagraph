CXX=g++ #-std=c++11
CXXFLAGS+= -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
#CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-seqan-v2.0.0/include 
CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include
LDFGLAS=-lrt

default: all

all: metagraph

metagraph: metagraph.o
	$(CXX) $(LDFLAGS) -o metagraph metagraph.o

metagraph.o: metagraph.cpp
	$(CXX) $(CXXFLAGS) -c -o metagraph.o $(LDFLAGS) metagraph.cpp

clean:
	rm -f metagraph.o metagraph

.PHONY: default all clean
