#include <map>
#include <string>
#include <vector>

#include "Models.hpp"

namespace dbg_database {

class MockDatabase //: public IDatabase
{

private:
    std::map<std::string, std::vector<std::string> > datastore;
    template<typename T> std::vector<T> concat(std::vector<T> vec1, std::vector<T> vec2);

public:
    MockDatabase();
    
    /**
     * Does all the necessary work to associate a given k-mer with an
     * annotation
     */
     int annotate_kmer(models::Kmer kmer, models::Annotation annotation);

     int annotate_kmer(std::string kmer, std::string raw_tag);

    /**
     * Looks up the annotation for a kmer.
     */
    models::Annotation get_annotation(models::Kmer kmer);

    models::Annotation get_annotation(std::string raw_kmer);

    /**
     * N.B. not in IDatabase.
     */
    std::map<std::string, std::vector<std::string> > getDatastore();
};

// MockDatabase::MockDatabase() {};
}
