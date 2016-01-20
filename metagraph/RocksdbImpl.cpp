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

    void partition(std::vector<std::string> columns) {
        rocksdb::Options options;
        options.create_if_missing = true;

        std::vector<rocksdb::ColumnFamilyDescriptor> column_families;
        column_families.push_back(rocksdb::ColumnFamilyDescriptor(rocksdb::kDefaultColumnFamilyName,
                                                                           rocksdb::ColumnFamilyOptions()));
        // make the partition
        for (auto it = columns.begin(); it < columns.end(); ++it) {
            rocksdb::ColumnFamilyHandle* cf;
            status = db->CreateColumnFamily(rocksdb::ColumnFamilyOptions(), *it, &cf);
            column_families.push_back(rocksdb::ColumnFamilyDescriptor(*it, rocksdb::ColumnFamilyOptions()));
            delete cf;
        }

        delete db;

        // reopen with new column families
        std::vector<rocksdb::ColumnFamilyHandle*> handles;
        status = rocksdb::DB::Open(options, dbpath, column_families, &handles, &db);
        assert(status.ok());
    }

public:
    RocksdbImpl(std::string the_dbpath) : dbpath(the_dbpath) {
        rocksdb::Options options;
        options.create_if_missing = true;
        status = rocksdb::DB::Open(options, dbpath, &db);

        std::vector<std::string> existing_column_families;
        rocksdb::DB::ListColumnFamilies(options, dbpath, &existing_column_families);
        // if (existing_column_families.size() == 1) {
        //     std::vector<std::string> column_names;
        //     column_names.push_back(ANNOTATION_COLUMN);
        //     column_names.push_back(KMER_COLUMN);
        //     partition(column_names);
        // }
    };

    void annotate_kmer(std::string kmer, std::string tag) {
        std::stringstream key;
        key << kmer << ":" << tag;

        std::string ignored = "";

        status = db->Put(rocksdb::WriteOptions(), key.str(), ignored);
        assert(status.ok());
    }

    std::string get_annotation(std::string kmer) {
        std::string ignored;
        db->Get(rocksdb::ReadOptions(), kmer, &ignored);

        std::cout << "get_annotation: " << kmer << std::endl;
        std::vector<std::string> ret;
        ret.push_back("FOOBAR");

        auto it = db->NewIterator(rocksdb::ReadOptions());
        it->Seek(kmer);
        std::cout << it->key().ToString() << std::endl;
        // iter->Next();
        // for (iter->Seek(kmer); iter->Valid(); iter->Next()) {
        // std::cout << "within iter" << std::endl;
        //     // do something
        //     iter->key();
        // }
        
        return ret[0];
    }
};

} // namespace dbg_database
