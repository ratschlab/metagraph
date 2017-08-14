LIBMAUS_DIR=external-libraries/libmaus2-lib
ROCKSDB_DIR=external-libraries/rocksdb
SDSL_DIR=external-libraries/sdsl

CXX=g++ -std=c++11 -fopenmp -D_THREAD_SAFE -pthread -O3 -g -msse2  #-march=core-avx-i -mtune=core-avx-i -mavx2 #-gdwarf-3 #-pg
CXXFLAGS+= -W -Wall -Wno-long-long -Wno-variadic-macros #-pedantic
CXXFLAGS+= -I$(LIBMAUS_DIR)/include -I. -L$(LIBMAUS_DIR)/lib -I$(SDSL_DIR)/include -L$(SDSL_DIR)/lib #-I$(ROCKSDB_DIR)/include -L$(ROCKSDB_DIR)
LDFLAGS=-lrt -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256 -lsdsl -lhts -lz #-lrocksdb -lsnappy -lbz2

### Makefile generated with ideas from https://stackoverflow.com/questions/2394609/makefile-header-dependencies

### all sources
SOURCES = $(wildcard *.cpp)

BINARY = metagraph
BUILD_DIR = ./build

### generate objects from all sources
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)
### generate dependencies from objects
DEPS = $(OBJECTS:%.o=%.d)

default: $(BINARY)

#$(BINARY): $(BUILD_DIR)/$(BINARY)

$(BINARY) : $(OBJECTS)
	### Create build directories with structure according to sources
	mkdir -p $(@D)
	### Link all object files
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

### include all dependency files
-include $(DEPS)

### build targets for single object files
$(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(@D)
	### -MMD creates a dependency file (.d) for each object file (.o)
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

bitest: bitest.cpp
	g++ -g -O3 -D__STDC_LIMIT_MACROS bitest.cpp -lhts -lz

clean:
	rm -rf $(BINARY) $(BUILD_DIR)

.PHONY: default all clean
