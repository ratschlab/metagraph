#include "kmer.hpp"


bool KMer::compare_kmer_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1.seq >> ((minus + 1) * kBitsPerChar)
             == k2.seq >> ((minus + 1) * kBitsPerChar);
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
    return ((seq >> (kBitsPerChar * i)) % kMax).convert_to<TAlphabet>();
}

std::ostream& operator<<(std::ostream &os, const KMer &kmer) {
    return os << kmer.seq;
}
