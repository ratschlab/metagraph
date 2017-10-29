#ifndef __DBG_SUCC_BOOST__
#define __DBG_SUCC_BOOST__


#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include "kseq.h"

#define BPC 3
typedef boost::multiprecision::uint256_t ui256;
typedef boost::multiprecision::cpp_int cpp_int;

namespace kmer_boost {

    typedef struct __KMer {
        ui256 seq;    // kmer sequences
       // cpp_int annot;  // annotation vector

        
        std::vector<uint32_t> annot;
        __KMer() {}
        __KMer(ui256 seq_) : seq(seq_) {}
        __KMer(ui256 seq_, size_t annot_):
            seq(seq_) {
                if (annot_ > 0)
                    annot.push_back(annot_ - 1);
            }
       
    } KMer;
    bool operator<(const KMer& a, const KMer& b);
    
    typedef struct __KMer_anno : public __KMer {
        ui256 seq;    // kmer sequences
        //cpp_int annot;  // annotation vector
        
        std::vector<uint32_t> annot;
        __KMer_anno(ui256 seq_) : seq(seq_) {}
        __KMer_anno(ui256 seq_, size_t annot_):
            seq(seq_) {
                if (annot_ > 0)
                    annot.push_back(annot_ - 1);
            }
        
    } KMer_anno;
    bool operator<(const KMer_anno& a, const KMer_anno& b);
    
    bool compare_kmer_suffix(const ui256 &k1, const ui256 &k2, const uint64_t &minus=0);
    
    ui256 stokmer(const char *seq, const uint64_t &len, const char *nt_lookup);
    
    ui256 stokmer(const kstring_t &seq, const char *nt_lookup);
    
    ui256 stokmer(const std::string &seq, const char *nt_lookup);
    
    uint8_t getW(ui256 kmer);
    
    //zero-based
    //if i==-1, then it gets W as a char
    char getPos(ui256 &kmer, int64_t i, const std::string &alphabet, const uint64_t &alph_size);
    
    char* kmertos(ui256 kmer, const std::string &alphabet, const uint64_t &alph_size);
    
}
#endif
