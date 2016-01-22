#include <assert.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <set>

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
    std::string ANNOTATION_COLUMN = "annotation";
    std::string KMER_COLUMN = "kmer";

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

        std::vector<std::string> existing_column_families;
        rocksdb::DB::ListColumnFamilies(options, dbpath, &existing_column_families);
    };

    void annotate_kmer(std::string kmer, std::string tag) {
        std::stringstream key;
        key << kmer << ":" << tag;

        std::string ignored = "";

        status = db->Put(rocksdb::WriteOptions(), key.str(), ignored);
        assert(status.ok());
    }

    std::vector<std::string> get_annotation(std::string kmer) {
        std::string ignored;
        db->Get(rocksdb::ReadOptions(), kmer, &ignored);

        auto it = db->NewIterator(rocksdb::ReadOptions());
        it->Seek(kmer);

        std::vector<std::string> ret;
        if (! it->Valid()) {
            return ret;
        }

        for (it->Seek(kmer);
             it->Valid() && (it->key().ToString().compare(0, kmer.length(), kmer) == 0);
             it->Next()) {
            ret.push_back(it->key().ToString().substr(kmer.length()+1));
        }
        
        return ret;
    }
};

} // namespace dbg_database
