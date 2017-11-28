#ifndef __DBG_SUCC_BOOST__
#define __DBG_SUCC_BOOST__

#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>

#include <htslib/kseq.h>

typedef boost::multiprecision::uint256_t ui256;
typedef boost::multiprecision::cpp_int cpp_int;

const int kBPC = 3;


class KMer {
    ui256 seq;    // kmer sequences

  public:
    std::vector<uint32_t> annot;
    KMer() {}
    explicit KMer(const ui256 &seq_) : seq(seq_) {}
    KMer(const KMer &k, size_t annot_) : KMer(k.seq, annot_) {}
    KMer(ui256 seq_, size_t annot_) : seq(seq_) {
            if (annot_ > 0)
                annot.push_back(annot_ - 1);
        }

    friend std::ostream& operator<<(std::ostream &os, const KMer &kmer) {
        os << kmer.seq;
        return os;
    }

    uint8_t getW() const;

    bool operator<(const KMer &other) const { return seq < other.seq; }
    bool operator==(const KMer &other) const { return seq == other.seq; }

    static bool compare_kmer_suffix(const KMer &k1, const KMer &k2, uint64_t minus = 0);

    template <class Map>
    static KMer stokmer(const char *seq, const uint64_t len, Map to_alphabet) {
        if (len * kBPC >= 256 || len < 2) {
            std::cerr << "String must be between lengths 2 and " << 256 / kBPC << ".\n";
            exit(1);
        }

        ui256 kmer = 0;
        uint8_t cur;
        uint8_t maxn = 1 << kBPC;

        for (int i = len - 2; i >= 0; --i) {
            kmer <<= kBPC;
            //this makes sure the range is 1-7 instead of 0-6
            cur = to_alphabet(seq[i]) + 1;
            if (cur < maxn) {
                kmer += cur;
            } else {
                std::cerr << "Alphabet size too big for the given number of bits. Alphabet size must be < " << maxn << "\n";
                exit(1);
            }
        }
        kmer <<= kBPC;
        kmer += to_alphabet(seq[len - 1]) + 1;
        return KMer {kmer};
    }

    template <class Map>
    static KMer stokmer(const std::string &seq, Map to_alphabet) {
        return stokmer(seq.c_str(), seq.length(), to_alphabet);
    }

    //zero-based
    //if i==-1, then it gets W as a char
    static char getPos(const ui256 &kmer, int64_t i, const std::string &alphabet, const uint64_t &alph_size);

    static char getPos(const KMer &kmer, int64_t i, const std::string &alphabet, const uint64_t &alph_size) {
        return getPos(kmer.seq, i, alphabet, alph_size);
    }

    static char* kmertos(const KMer &k, const std::string &alphabet, const uint64_t &alph_size);
};

#endif // __DBG_SUCC_BOOST__
