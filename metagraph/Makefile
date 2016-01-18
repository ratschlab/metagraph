LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB=/home/gideon/Scratch/rocksdb
SEQAN_DIR=external-libraries/seqan-include

#TODO revert these changes.
CXX=g++ -std=c++11 -fopenmp #-O3 #-pg
#TODO bring back Makefile to previous highly sensitive state.
#CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(SEQAN_DIR) -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib
LDFLAGS=-lrt -lmaus2

default: all

all: metagraph

debug: metagraph_debug

# TODO cleanup rocks db deps.
metagraph: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -I $(ROCKSDB)/include -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp MockDatabase.hpp MockDatabase.cpp IDatabaseImpl.cpp $(ROCKSDB)/librocksdb.a -lpthread -lrt -lsnappy -lgflags -lz -lbz2 -llz4 -ljemalloc -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native -std=c++11  -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp
	$(CXX) $(CXXFLAGS) -o metagraph metagraph.cpp dbg_succinct_libmaus.cpp $(LDFLAGS) -D DBGDEBUG

# TODO database task that creates a database.o object file.

database:
	$(CXX) $(CXXFLAGS) -I $(ROCKSDB)/include IDatabaseImpl.cpp DatabaseCli.cpp -o database-cli $(ROCKSDB)/librocksdb.a -lpthread -lrt -lsnappy -lgflags -lz -lbz2 -llz4 -ljemalloc -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native -std=c++11  -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native $(LDFLAGS)

database_test:
	$(CXX) $(CXXFLAGS) TestDatabase.cpp MockDatabase.hpp MockDatabase.cpp $(LDFLAGS)
	./a.out
	rm -f a.out

clean:
	rm -f a.out
	rm -f metagraph.o metagraph
	rm -f IDatabase.hpp.* IDatabaseImpl.cpp.* IDatabaseImpl.o Models.hpp.*
	rm -f *.gch
	rm -f database-cli

.PHONY: default all clean
