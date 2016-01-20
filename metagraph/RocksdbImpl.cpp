#include <assert.h>
#include <iostream>
#include <stdexcept>
#include <string>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

// TODO what is the difference between `#include ".*"` and `#include <.*>"`
#include "IDatabase.hpp"

namespace dbg_database {

class RocksdbImpl : public IDatabase {

private:
    std::string dbpath;
    rocksdb::DB* db;
    rocksdb::Status status;

    template<typename T> std::vector<T> concat(std::vector<T> vec1, std::vector<T> vec2) {
        // http://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
        std::vector<T> ret;
        ret.reserve( vec1.size() + vec2.size() ); // preallocate memory
        ret.insert( ret.end(), vec1.begin(), vec1.end() );
        ret.insert( ret.end(), vec2.begin(), vec2.end() );

        return ret;
    }

public:
    RocksdbImpl(std::string the_dbpath) : dbpath(the_dbpath) {
        rocksdb::Options options;
        options.create_if_missing = true;
        status = rocksdb::DB::Open(options, dbpath, &db);
    };

    void annotate_kmer(std::string kmer, std::string tag) {
        assert(status.ok());
        db->Put(rocksdb::WriteOptions(), kmer, tag);
        assert(status.ok());
    }

    std::string get_annotation(std::string kmer) {
        std::string ret;
        assert(status.ok());
        db->Get(rocksdb::ReadOptions(), kmer, &ret);
        assert(status.ok());
        return ret;
    }
};

} // namespace dbg_database
