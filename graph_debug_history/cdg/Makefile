CXX=gcc -g # -O2
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros -fpermissive
#CXXFLAGS+= -I/cbio/grlab/home/akahles/seqan/seqan-1.4.2/core/include -I /cbio/grlab/home/akahles/software/cpp-btree-1.0.1 -I/cbio/grlab/home/akahles/software/sdsl-lite-lib/include -I. -L/cbio/grlab/home/akahles/software/sdsl-lite-lib/lib
LDFLAGS= #-lrt -lsdsl -ldivsufsort -ldivsufsort64

default: all

all: sdg

sdg: sdg.c
	$(CXX) $(CXXFLAGS) -o sdg sdg.c $(LDFLAGS)

clean:
	rm -f sdg

.PHONY: default all clean
