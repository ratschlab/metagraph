LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB_DIR=external-libraries/rocksdb

CXX=g++ -std=c++11 -fopenmp -D_THREAD_SAFE -std=c++0x -pthread -O3 #-g -gdwarf-3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib -I$(ROCKSDB_DIR)/include -L$(ROCKSDB_DIR)
LDFLAGS=-lrt -lmaus2 -lrocksdb -lsnappy -lbz2

default: all

all: metagraph

debug: metagraph_debug

metagraph: dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp metagraph.cpp config.cpp serialization.hpp datatypes.hpp
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS)

metagraph_debug: metagraph.cpp dbg_succinct_libmaus.cpp dbg_succinct_libmaus.hpp config.hpp serialization.hpp datatypes.hpp 
	$(CXX) $(CXXFLAGS) -o metagraph dbg_succinct_libmaus.cpp metagraph.cpp config.cpp $(LDFLAGS) -D DBGDEBUG

# TODO database task that creates a database.o object file.

database:
	$(CXX) $(CXXFLAGS) -I $(ROCKSDB)/include RocksdbImpl.cpp DatabaseCli.cpp -o database-cli $(ROCKSDB)/librocksdb.a -lpthread -lrt -lsnappy -lgflags -lz -lbz2 -llz4 -ljemalloc -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native -std=c++11  -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native $(LDFLAGS)

database_test: 
	$(CXX) $(CXXFLAGS) -I $(ROCKSDB)/include TestDatabase.cpp MockDatabase.cpp RocksdbImpl.cpp dbg_succinct_libmaus.cpp $(ROCKSDB)/librocksdb.a -lpthread -lrt -lsnappy -lgflags -lz -lbz2 -llz4 -ljemalloc -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native -std=c++11  -DROCKSDB_PLATFORM_POSIX -DROCKSDB_LIB_IO_POSIX  -DOS_LINUX -fno-builtin-memcmp -DROCKSDB_FALLOCATE_PRESENT -DSNAPPY -DGFLAGS=gflags -DZLIB -DBZIP2 -DLZ4 -DROCKSDB_MALLOC_USABLE_SIZE -march=native $(LDFLAGS) 
	./a.out
	rm -f a.out

clean:
	rm -f a.out
	rm -f metagraph.o metagraph
	rm -f IDatabase.hpp.* RocksdbImpl.cpp.* RocksdbImpl.o Models.hpp.*
	rm -f *.gch
	rm -f database-cli

.PHONY: default all clean
