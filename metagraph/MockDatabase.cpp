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

#include <IDatabase.hpp>

// TODO this type definition is polluting everything.
typedef seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'X'> > Dna5F; 

namespace dbg_database {

class MockDatabase : public IDatabase {

private:
    std::map<std::string, std::string> datastore;

public:
    MockDatabase() {};
    
    // TODO handle multiple tag resolution: dedup and sort.
    void annotate_kmer(std::string kmer, std::string tag) {
        datastore[kmer] = tag;
    };

    std::vector<std::string> get_annotation(std::string kmer) {
        std::vector<std::string> ret;
        ret.push_back(kmer);
        return ret;
    };

    /**
     * N.B. not in IDatabase.
     */
    std::map<std::string, std::string> getDatastore() {
        return datastore;
    }
};

} // namespace dbg_database
