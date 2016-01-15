LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB=/home/gideon/Scratch/rocksdb
SEQAN_DIR=external-libraries/seqan-include

#TODO revert these changes.
CXX=g++ -std=c++11 -fopenmp #-O3 #-pg
#TODO bring back Makefile to previous highly sensitive state.
#CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(SEQAN_DIR) -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib
LDFLAGS=-lrt -lmaus2

ROCKSDB_LDFLAGS= -lpthread -lrt -lsnappy -lgflags -lz -lbz2
ROCKSDB_PLATFORM_CXXFLAGS=-std=c++11  -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=google -DZLIB -DBZIP2 -DROCKSDB_MALLOC_USABLE_SIZE -march=native 

default: all

all: metagraph

debug: metagraph_debug

metagraph: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp MockDatabase.hpp MockDatabase.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp $(LDFLAGS) -D DBGDEBUG

# TODO database task that creates a database.o object file.

database:
	$(CXX) $(CXXFLAGS) -I $(ROCKSDB)/include IDatabaseImpl.cpp $(ROCKSDB)/librocksdb.a $(ROCKSDB_LDFLAGS) $(ROCKSDB_PLATFORM_CXXFLAGS) $(LDFLAGS)

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
