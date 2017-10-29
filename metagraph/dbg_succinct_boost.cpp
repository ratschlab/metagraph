#include "dbg_succinct_boost.hpp"


bool kmer_boost::operator<(const kmer_boost::KMer& a, const kmer_boost::KMer& b) {
    return a.seq < b.seq;
}

bool kmer_boost::operator<(const kmer_boost::KMer_anno& a, const kmer_boost::KMer_anno& b) {
    if (a.seq == b.seq) {
        //return a.annot < b.annot;
        return a.annot.size() < b.annot.size();
    } else {
        return a.seq < b.seq;
    }
}

bool kmer_boost::compare_kmer_suffix(const ui256 &k1, const ui256 &k2, const uint64_t &minus) {
    return (k1 >> ((minus+1) * BPC)) == (k2 >> ((minus+1)* BPC));
}

ui256 kmer_boost::stokmer(const char *seq, const uint64_t &len, const char *nt_lookup) {
    if (len * BPC >= 256 || len < 2) {
        std::cerr << "String must be between lengths 2 and " << 256/BPC << ".\n";
        exit(1);
    }   
    ui256 kmer = 0;
    uint8_t cur,maxn=(1<<BPC);
    for (int i = len-2; i>=0;--i) {
        kmer <<= BPC;
        //this makes sure the range is 1-7 instead of 0-6
        cur = nt_lookup[(uint8_t)seq[i]]+1;
        if (cur < maxn) {
            kmer += cur;
        } else {
            std::cerr << "Alphabet size too big for the given number of bits. Alphabet size must be < " << maxn << "\n";
            exit(1);
        }
    }   
    kmer <<= BPC;
    kmer += nt_lookup[(uint8_t)seq[len-1]]+1;
    return kmer;
}

ui256 kmer_boost::stokmer(const kstring_t &seq, const char *nt_lookup) {
    return kmer_boost::stokmer(seq.s, seq.l, nt_lookup);
}

ui256 kmer_boost::stokmer(const std::string &seq, const char *nt_lookup) {
    return kmer_boost::stokmer(seq.c_str(), seq.length(), nt_lookup);
}

uint8_t kmer_boost::getW(ui256 kmer) {
    uint8_t curc=0;
    for (size_t i=0;i<BPC;i++) {
        curc += (kmer % 2).convert_to<uint8_t>() << i;
        kmer >>= 1;
    }
    if (!curc)
        return 127;
    return curc-1;
}

//zero-based
//if i==-1, then it gets W as a char
char kmer_boost::getPos(ui256 &kmer, int64_t i, const std::string &alphabet, const uint64_t &alph_size) {
    uint8_t curw = kmer_boost::getW(kmer >> ((BPC)*(i+1)));
    if (curw == 127)
        return 0;
    return (curw < alph_size ? alphabet[curw] : 'N');
}

char* kmer_boost::kmertos(ui256 kmer, const std::string &alphabet, const uint64_t &alph_size) {
    char *seq = (char*)calloc(256/BPC+1, sizeof(char));
    for (uint8_t curpos=0;kmer;++curpos, kmer >>= BPC) {
        seq[curpos] = kmer_boost::getPos(kmer, -1, alphabet, alph_size);
    }
    return seq;
}

