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

using namespace std;
using namespace rocksdb;

class IDatabaseImpl
{
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

    // virtual int annotate_kmer(models::Kmer kmer, models::Annotation annotation) = 0;
    // virtual models::Annotation get_annotation(models::Kmer kmer) = 0;
};

int main(int argc, char const ** argv) {
    IDatabaseImpl database;

    database.annotate_kmer("TCGA", "cancer");
    cout << database.get_annotation("TCGA") << endl;

    // // open DB
    // rocksdb::Options options;
    // options.create_if_missing = true;

    // DB* db;
    // Status status = DB::Open(options, kDBPath, &db);
    // assert(status.ok());

    // std::string value;

    // // db->Put(WriteOptions(), "akey", value);
    // db->Get(ReadOptions(), "akey", &value);

    // cout << value << endl;
    
    // // rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key1, &value);
    // // if (s.ok()) s = db->Put(rocksdb::WriteOptions(), key2, value);
    // // if (s.ok()) s = db->Delete(rocksdb::WriteOptions(), key1);

    // // close DB
    // delete db;
}

// class IDatabaseImpl //: public IDatabase
// {

// private:

//     /**
//      * Adds tags that we don't know about yet. Ignore ones that we do.
//      */
//     void update_tags(pqxx::connection* conn, Annotation annotation) {
//         for (size_t i=0; i < annotation.tags.size(); i++) {
//             Tag tag = annotation.tags[i];
//             try {
//                 pqxx::work transaction(*conn, "annotate_kmer");

//                 std::string query = "INSERT INTO tags (tag) VALUES ('";
//                 query += tag.name;
//                 query += "')";
//                 transaction.exec(query);
//                 transaction.commit();
//             } catch (pqxx::unique_violation ignored) {
//             }
//         }
//     };

//     std::vector<int> get_tag_ids(std::vector<Tag> tags) {
//         // "SELECT id FROM tags WHERE tag IN ()";
//     }

//     int get_annotation_id(pqxx::connection* conn, Annotation Annotation) {
//         return 42;
//     };

//     std::string join(std::vector<std::string> strings, char delim) {
//     }

//     // TODO this code is duplicated.
//     template<typename T> std::vector<T> concat(std::vector<T> vec1, std::vector<T> vec2) {
//         // http://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
//         std::vector<T> ret;
//         ret.reserve( vec1.size() + vec2.size() ); // preallocate memory
//         ret.insert( ret.end(), vec1.begin(), vec1.end() );
//         ret.insert( ret.end(), vec2.begin(), vec2.end() );

//         return ret;
//     }
    
// public:

//     IDatabaseImpl();

//     virtual Annotation get_annotation(Kmer kmer) {
//         seqan::Shape<seqan::Dna5, seqan::SimpleShape> shape;
//         seqan::resize(shape, kmer.name.size());
//         int thehash = seqan::hash(shape, begin(kmer.name));
//         return get_annotation(thehash);
//     }

//     virtual Annotation get_annotation(int hash) {
//         pqxx::connection conn(connect_string);

//         pqxx::work indicator(conn, "contains_hash");
//         pqxx::result has_kmer = indicator.exec("SELECT * FROM kmers WHERE hash = " + std::to_string(hash));
        
//         indicator.commit();
//         if (has_kmer.size() < 1) {
//             // TODO custom exception?
//             throw std::out_of_range("kmer not found in db.");
//         }

//         // TODO this is a JOIN.
//         std::string query = "SELECT tags.tag FROM tags WHERE tags.id \
//         IN (SELECT unnest(tags) FROM annotations WHERE annotations.id \
//         IN (SELECT annotation_id FROM kmers WHERE hash = " +
//         std::to_string(hash) + "))";

//         pqxx::work transaction(conn, "get_annotation");

//         pqxx::result result = transaction.exec(query);
//         transaction.commit();

//         std::vector<Tag> tags;
//         for (pqxx::result::size_type i = 0; i < result.size(); ++i) {
//             Tag tag;
//             tag.name = result[i][0].c_str();
//             tags.push_back(tag);
//         }

//         Annotation ret;
//         ret.tags = tags;
//         return ret;
//     }
    
//     virtual int annotate_kmer(Kmer kmer, Annotation annotation) {

//         Annotation stored_annotation;
//         try {
//             stored_annotation = get_annotation(kmer);
//         } catch (std::out_of_range &ignored) {
//             std::vector<Tag> tags;
//             stored_annotation.tags = tags;
//         }

//         auto merged = concat(annotation.tags, stored_annotation.tags);

//         std::cout << merged.size() << std::endl;

//         // TODO
//         // for (int i=0; i < merged.size(); ++i)
//         //     std::cout << merged[i] << std::endl;

//         // look up kmer
//         // update annotations
//         // - add new tags
//         // - check what annotations are missing
//         // - add new ones if necessary

//         // pqxx::connection conn(connect_string);
//         // update_tags(&conn, annotation);
//         // int annotation_id = get_annotation_id(&conn, annotation);
        
//         return 42;
//     };
// };

// IDatabaseImpl::IDatabaseImpl() { }
