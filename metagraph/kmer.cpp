#include "kmer.hpp"

KMer operator>>(const KMer &kmer, const size_t &shift) {
    KMer that(kmer);
    that >>= shift;
    return that;
}

bool KMer::compare_kmer_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1 >> ((minus + 1) * kBitsPerChar)
             == k2 >> ((minus + 1) * kBitsPerChar);
}

uint64_t KMer::operator[](size_t i) const {
    assert(get_digit<kBitsPerChar>(i) > 0);
    return get_digit<kBitsPerChar>(i) - 1;
}

std::string KMer::to_string(const std::string &alphabet) const {
    std::string seq;
    seq.reserve(256 / kBitsPerChar + 1);

    uint64_t cur;
    for (size_t i = 0; (cur = get_digit<kBitsPerChar>(i)); ++i) {
        seq.push_back(alphabet.at(cur - 1));
    }
    return seq;
}

std::ostream& operator<<(std::ostream &os, const KMer &kmer) {
    os << kmer.s_[0];
    for (uint8_t i = 1; i < 4; ++i) {
        os << " " << kmer.s_[i];
    }
    return os;
}

KMer operator<<(const KMer &kmer, const size_t &shift) {
    KMer that(kmer);
    that <<= shift;
    return that;
}
