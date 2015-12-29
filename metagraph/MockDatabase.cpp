#include <algorithm>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "MockDatabase.hpp"

using namespace seqan;
using namespace std;

// TODO this type definition is polluting everything.
typedef seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'X'> > Dna5F; 

namespace dbg_database {

// TODO this code is duplicated. Should probably go into some sort
// of utils file.
template<typename T> std::vector<T> MockDatabase::concat(std::vector<T> vec1, std::vector<T> vec2) {
    // http://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
    std::vector<T> ret;
    ret.reserve( vec1.size() + vec2.size() ); // preallocate memory
    ret.insert( ret.end(), vec1.begin(), vec1.end() );
    ret.insert( ret.end(), vec2.begin(), vec2.end() );

    return ret;
}

MockDatabase::MockDatabase() {};
    
/**
 * Does all the necessary work to associate a given k-mer with an
 * annotation
 */
int MockDatabase::annotate_kmer(models::Kmer kmer, models::Annotation annotation) {
    vector<string> existing_tags = datastore[kmer.name];

    // extract tags as strings
    vector<string> raw_tags;
    for (auto it = annotation.tags.begin(); it != annotation.tags.end(); ++it) {
        raw_tags.push_back(it->name);
    }

    vector<string> tags = concat(existing_tags, raw_tags);

    std::set<string> deduped(tags.begin(), tags.end());
    tags.assign(deduped.begin(), deduped.end());

    datastore[kmer.name] = tags;
    return datastore.size();
};

int MockDatabase::annotate_kmer(string kmer, string raw_tag) {
    vector<models::Tag> tags;
    tags.push_back((models::Tag) { raw_tag });
    return annotate_kmer((models::Kmer) {kmer}, (models::Annotation) {tags});
};

/**
 * Looks up the annotation for a kmer.
 */
models::Annotation MockDatabase::get_annotation(models::Kmer kmer) {
    vector<string> raw_tags = datastore[kmer.name];

    vector<models::Tag> tags;
    for (auto it = raw_tags.begin(); it != raw_tags.end(); ++it) {
        models::Tag tag;
        tag.name = *it;
        tags.push_back(tag);
    }

    models::Annotation to_return;
    to_return.tags = tags;
    return to_return;
};

     models::Annotation MockDatabase::get_annotation(string raw_kmer) {
        return get_annotation( (models::Kmer) { raw_kmer });
    }

    /**
     * N.B. not in IDatabase.
     */
    map<string, vector<string> > MockDatabase::getDatastore() {
        return datastore;
    }

} // namespace dbg_database
