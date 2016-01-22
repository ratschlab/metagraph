#ifndef __IDATABASE_HPP__
#define __IDATABASE_HPP__

#include <iostream>

namespace dbg_database {

class IDatabase
{
public:

    // TODO destructor??
    // // always inclue destructor in interface to help runtime do
    // // cleanup on the correct object.
    // virtual ~IDatabase();

    /**
     * Does all the necessary work to associate a given k-mer with an
     * annotation
     */
    virtual void annotate_kmer(std::string kmer, std::string raw_tag) = 0;

    /**
     * Looks up the annotation for a kmer.
     */
    virtual std::vector<std::string> get_annotation(std::string raw_kmer) = 0;
};

} // dbg_database

#endif
