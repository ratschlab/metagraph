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
    KMer(KMer &&other) : seq_(other.seq_) {}
    KMer(const KMer &other) : seq_(other.seq_) {}
    explicit KMer(ui256 &&seq) : seq_(seq) {}
    explicit KMer(const ui256 &seq) : seq_(seq) {}

    KMer& operator=(KMer &&other) { seq_ = other.seq_; return *this; }
    KMer& operator=(const KMer &other) { seq_ = other.seq_; return *this; }

    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    TAlphabet operator[](size_t i) const { return get(i) - 1; }

    static bool compare_kmer_suffix(const KMer &k1,
                                    const KMer &k2, size_t minus = 0);

    template <class Map, class String>
    static KMer from_string(const String &seq, Map &&to_alphabet) {
        assert(seq.size() * kBitsPerChar < 256 && seq.size() >= 2
                && "String must be between lengths 2 and 256 / kBitsPerChar");

        ui256 kmer = 0;
        for (int i = seq.size() - 2; i >= 0; --i) {
            uint8_t cur = to_alphabet(seq[i]) + 1;

            assert(cur < kMax && "Alphabet size too big for the given number of bits");

            kmer = (kmer << kBitsPerChar) + cur;
        }
        kmer <<= kBitsPerChar;
        kmer += to_alphabet(seq[seq.size() - 1]) + 1;
        return KMer(kmer);
    }

    std::string to_string(const std::string &alphabet) const;

  private:
    ui256 seq_; // kmer sequence

    TAlphabet get(size_t i) const;
};

#endif // __KMER_HPP__
