#include <map>
#include <string>
#include <vector>

namespace dbg_database {

class MockDatabase : public IDatabase
{

private:
    std::map<std::string, std::vector<std::string> > datastore;
    template<typename T> std::vector<T> concat(std::vector<T> vec1, std::vector<T> vec2);

public:
    MockDatabase();
    
    void annotate_kmer(std::string kmer, std::string raw_tag);

    string get_annotation(std::string raw_kmer);

    /**
     * N.B. not in IDatabase.
     */
    std::map<std::string, std::vector<std::string> > getDatastore();
};

// MockDatabase::MockDatabase() {};
} // namespace dbg_database
