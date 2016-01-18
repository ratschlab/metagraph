#include <assert.h>
#include <iostream>
#include <stdexcept>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

// TODO what is the difference between `#include ".*"` and `#include <.*>"`
#include "IDatabase.hpp"

using namespace std;
using namespace rocksdb;

class IDatabaseImpl : public IDatabase {

private:
    string kDBPath = "/tmp/rocksdb-test";
    DB* db;
    Status status;
    
public:
    IDatabaseImpl() {
        rocksdb::Options options;
        options.create_if_missing = true;
        status = DB::Open(options, kDBPath, &db);
    };

    void annotate_kmer(string raw_kmer, string raw_tag) {
        assert(status.ok());
        db->Put(WriteOptions(), raw_kmer, raw_tag);
        assert(status.ok());
    }

    string get_annotation(string raw_kmer) {
        string ret;
        assert(status.ok());
        db->Get(ReadOptions(), raw_kmer, &ret);
        assert(status.ok());
        return ret;
    }
};
