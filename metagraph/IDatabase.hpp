#include "Models.hpp"

// N.B. How do you do interfaces in C++? When you learn to do that, incorporate this interface.
class IDatabase
{
public:

    // always inclue destructor in interface to help runtime do
    // cleanup on the correct object.
    virtual ~IDatabase();

    /**
     * Does all the necessary work to associate a given k-mer with an
     * annotation
     */
    virtual int annotate_kmer(models::Kmer kmer, models::Annotation annotation) = 0;
    virtual int annotate_kmer(std::string kmer, std::string raw_tag) = 0;

    /**
     * Looks up the annotation for a kmer.
     */
    virtual models::Annotation get_annotation(models::Kmer kmer) = 0;
    virtual models::Annotation get_annotation(std::string raw_kmer) = 0;
};
