#include "kmer.hpp"


uint8_t KMer::getW() const {
    ui256 kmer = seq;
    uint8_t curc = 0;
    for (size_t i = 0; i < kBPC; ++i) {
        curc += (kmer % 2).convert_to<uint8_t>() << i;
        kmer >>= 1;
    }
    if (!curc)
        return 127;
    return curc - 1;
}

bool KMer::compare_kmer_suffix(const KMer &k1, const KMer &k2, uint64_t minus) {
    return k1.seq >> ((minus + 1) * kBPC) == k2.seq >> ((minus + 1) * kBPC);
}

//zero-based
//if i==-1, then it gets W as a char
char KMer::getPos(const ui256 &kmer, int64_t i, const std::string &alphabet, const uint64_t &alph_size) {
    uint8_t curw = KMer(kmer >> (kBPC * (i + 1))).getW();
    if (curw == 127)
        return 0;
    return (curw < alph_size ? alphabet[curw] : 'N');
}

char* KMer::kmertos(const KMer &k, const std::string &alphabet, const uint64_t &alph_size) {
    ui256 kmer = k.seq;
    char *seq = (char*)calloc(256 / kBPC + 1, sizeof(char));
    for (uint8_t curpos = 0; kmer; ++curpos, kmer >>= kBPC) {
        seq[curpos] = getPos(kmer, -1, alphabet, alph_size);
    }
    return seq;
}

