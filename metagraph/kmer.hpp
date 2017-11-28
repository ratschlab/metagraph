#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <algorithm>

#include <htslib/kseq.h>

typedef boost::multiprecision::uint256_t ui256;
typedef boost::multiprecision::cpp_int cpp_int;
typedef uint64_t TAlphabet;

const int kBitsPerChar = 3;
const int kMax = 1llu << kBitsPerChar;


class KMer {
    friend std::ostream& operator<<(std::ostream &os, const KMer &kmer);

  public:
    std::vector<uint32_t> annot;

    KMer() {}

    explicit KMer(const ui256 &seq_) : seq(seq_) {}

    KMer(const KMer &k, size_t annot_) : KMer(k.seq, annot_) {}

    KMer(ui256 seq_, size_t annot_) : seq(seq_) {
        if (annot_ > 0)
            annot.push_back(annot_ - 1);
    }

    TAlphabet operator[](size_t i) const { return get(i) - 1; }

    bool operator<(const KMer &other) const { return seq < other.seq; }

    bool operator==(const KMer &other) const { return seq == other.seq; }

    static bool compare_kmer_suffix(const KMer &k1,
                                    const KMer &k2, size_t minus = 0);

    template <class Map>
    static KMer from_string(const std::string &seq, Map &&to_alphabet) {
        assert(seq.length() * kBitsPerChar < 256 && seq.length() >= 2
                && "String must be between lengths 2 and 256 / kBitsPerChar");

        ui256 kmer = 0;
        for (int i = seq.length() - 2; i >= 0; --i) {
            uint8_t cur = to_alphabet(seq[i]) + 1;

            assert(cur < kMax && "Alphabet size too big for the given number of bits");

            kmer = (kmer << kBitsPerChar) + cur;
        }
        kmer <<= kBitsPerChar;
        kmer += to_alphabet(seq[seq.length() - 1]) + 1;
        return KMer(kmer);
    }

    std::string to_string(const std::string &alphabet) const;

  private:
    ui256 seq; // kmer sequence

    TAlphabet get(size_t i) const;
};

#endif // __KMER_HPP__
