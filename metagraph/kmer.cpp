#include "kmer.hpp"

KMer operator>>(const KMer &kmer, const size_t &shift) {
    KMer that(kmer);
    that >>= shift;
    return that;
}

size_t operator|(const KMer &kmer, const int &a) {
    KMer that(kmer);
    that |= a;
    return that.nth(0);
}


bool KMer::compare_kmer_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1 >> ((minus + 1) * kBitsPerChar)
             == k2 >> ((minus + 1) * kBitsPerChar);
}

std::string KMer::to_string(const std::string &alphabet) const {
    std::string seq;
    seq.reserve(256 / kBitsPerChar + 1);

    TAlphabet cur;
    for (size_t i = 0; (cur = get(i)); ++i) {
        seq.push_back(alphabet.at(cur - 1));
    }
    return seq;
}

TAlphabet KMer::get(size_t i) const {
    int bit = kBitsPerChar * i;
    if (bit % 64 == 0) {
        return this->nth(bit / 64) % kMax;
    } else {
        TAlphabet result = this->nth(bit / 64) >> (bit % 64);
        if (bit / 64 < 3) {
            return (result | ((this->nth(bit / 64 + 1)) << (64 - (bit % 64)))) % kMax;
        } else {
            return result % kMax;
        }
    }
}

std::ostream& operator<<(std::ostream &os, const KMer &kmer) {
    for (uint8_t i = 0; i < 4; ++i) {
        os << kmer.s_[i] << " ";
    }
    return os;
}

KMer operator<<(const KMer &kmer, const size_t &shift) {
    KMer that(kmer);
    that <<= shift;
    return that;
}

