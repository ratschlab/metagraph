SEQAN_DIR=external-libraries/seqan-include
LIBMAUS_DIR=external-libraries/libmaus2-lib

CXX=g++ -std=c++11 -fopenmp #-O3 #-pg
#TODO bring back Makefile to previous highly sensitive state.
#CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(SEQAN_DIR) -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib
LDFLAGS=-lrt -lmaus2

default: all

all: metagraph

debug: metagraph_debug

metagraph: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp MockDatabase.hpp MockDatabase.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp $(LDFLAGS) -D DBGDEBUG

# TODO database task that creates a database.o object file.

database:
	$(CXX) $(CXXFLAGS) IDatabaseImpl.cpp $(LDFLAGS)

database_test:
	$(CXX) $(CXXFLAGS) TestDatabase.cpp MockDatabase.hpp MockDatabase.cpp $(LDFLAGS)
	./a.out
	rm -f a.out

clean:
	rm -f a.out
	rm -f metagraph.o metagraph
	rm -f IDatabase.hpp.* IDatabaseImpl.cpp.* IDatabaseImpl.o Models.hpp.*
	rm -f *.gch

.PHONY: default all clean
