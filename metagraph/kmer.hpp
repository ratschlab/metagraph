#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <algorithm>
#include <vector>

#include <htslib/kseq.h>

typedef boost::multiprecision::uint256_t ui256;
typedef boost::multiprecision::cpp_int cpp_int;
typedef uint64_t TAlphabet;

const int kBitsPerChar = 3;
const int kMax = 1llu << kBitsPerChar;


class KMer {
    friend std::ostream& operator<<(std::ostream &os, const KMer &kmer);

  public:
    template <typename T>
    KMer(const T &arr, size_t k);

    template <class Map, class String>
    KMer(const String &seq, Map &&to_alphabet);

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

    std::string to_string(const std::string &alphabet) const;

  private:
    ui256 seq_; // kmer sequence

    TAlphabet get(size_t i) const;
};

template <typename T>
KMer::KMer(const T &arr, size_t k) : seq_(0) {
    assert(k * kBitsPerChar < 256 && k >= 2
            && "String must be between lengths 2 and 256 / kBitsPerChar");

    for (int i = k - 2; i >= 0; --i) {
        assert(arr[i] + 1 < kMax && "Alphabet size too big for the given number of bits");

        seq_ <<= kBitsPerChar;
        seq_ += arr[i] + 1;
    }
    seq_ <<= kBitsPerChar;
    seq_ += arr[k - 1] + 1;
}

template <class Map, class String>
KMer::KMer(const String &seq, Map &&to_alphabet) {
    std::vector<uint8_t> arr(seq.size());
    std::transform(seq.begin(), seq.end(), arr.begin(), to_alphabet);
    *this = KMer(arr.data(), arr.size());
}

#endif // __KMER_HPP__
