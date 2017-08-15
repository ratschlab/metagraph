#ifndef __DBG_SUCC_BOOST__
#define __DBG_SUCC_BOOST__


#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <vector>
#include <iterator>
#include "kseq.h"

#define BPC 3
typedef boost::multiprecision::uint256_t ui256;

bool compare_kmer_suffix(const ui256 &k1, const ui256 &k2, const uint64_t &minus=0) {
    return (k1 >> ((minus+1) * BPC)) == (k2 >> ((minus+1)* BPC));
}

ui256 stokmer(const char *seq, const uint64_t &len, const char *nt_lookup) {
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

ui256 stokmer(const kstring_t &seq, const char *nt_lookup) {
    return stokmer(seq.s, seq.l, nt_lookup);
}

ui256 stokmer(const std::string &seq, const char *nt_lookup) {
    return stokmer(seq.c_str(), seq.length(), nt_lookup);
}

uint8_t getW(ui256 kmer) {
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
char getPos(ui256 &kmer, int64_t i, const char *alphabet, const uint64_t &alph_size) {
    uint8_t curw = getW(kmer >> ((BPC)*(i+1)));
    if (curw == 127)
        return 0;
    return (curw < alph_size ? alphabet[curw] : 'N');
}

char* kmertos(ui256 kmer, const char *alphabet, const uint64_t &alph_size) {
    char *seq = (char*)calloc(256/BPC+1, sizeof(char));
    for (uint8_t curpos=0;kmer;++curpos, kmer >>= BPC) {
        seq[curpos] = getPos(kmer, -1, alphabet, alph_size);
    } 
    return seq;
}

size_t seqtokmer(std::vector<ui256> &kmers, const char *seq, const uint64_t &len, const uint64_t &k, const char *nt_lookup, bool add_bridge=true, unsigned int parallel=1, std::string suffix="") {
    char *bridge = (char*)malloc(k+2);
    memset(bridge, 'X', k); 
    bridge[k] = seq[0];
    bridge[k+1] = 0;
    size_t i=0;
    size_t nkmers=0;
    //size_t j=kmers.size();
    //std::cout << "Loading next sequence with " << parallel << " threads\n";
    //kmers.resize(kmers.size()+len-k+(2*k*add_bridge));
    //std::cerr << "Memory reserved\n";
    if (add_bridge) {
        for (i=0;i<std::min(k,len); ++i) {
            //kmers[j++]=stokmer(bridge, k+1, nt_lookup);
            if (std::string(bridge+k-suffix.length(),bridge+k) == suffix) {
                kmers.push_back(stokmer(bridge, k+1, nt_lookup));
            }
            nkmers++;
            memmove(bridge, bridge+1, k); 
            bridge[k]= (i+1 < len) ? seq[i+1] : 'X';
        }
    }
    if (k<len) {
        #pragma omp parallel num_threads(parallel)
        {
            std::vector<ui256> kmer_priv;
            #pragma omp for nowait
            for (i=0; i < len-k; ++i) {
                //if (!((i+1) % 1000)) {
                //    std::cout << "." << std::flush;
                    //if (!((i+1) % 100000)) {
                    //    std::cout << i+1 << "\n";
                    //} 
                //}
                //kmers[i+j]=stokmer(seq+i, k+1, nt_lookup);
                if (std::string(seq+i+k-suffix.length(), seq+i+k) == suffix) {
                    //std::cerr << std::string(seq+i, seq+i+k+1) << "\n";
                    kmer_priv.push_back(stokmer(seq+i, k+1, nt_lookup));
                }
            }
            #pragma omp critical
            kmers.insert(kmers.end(), std::make_move_iterator(kmer_priv.begin()), std::make_move_iterator(kmer_priv.end()));
        }
        nkmers += len-k;
        memcpy(bridge, seq+len-k, k); 
        bridge[k]='X';
        //j += len-k;
    }   
    if (add_bridge) {
        for (i=0;i<k;++i) {
            //kmers[j++]=stokmer(bridge, k+1, nt_lookup);
            if (std::string(bridge+k-suffix.length(),bridge+k) == suffix) {
                //std::cerr << std::string(bridge, bridge+k+1) << "\n";
                kmers.push_back(stokmer(bridge, k+1, nt_lookup));
            }
            nkmers++;
            memmove(bridge, bridge+1, k); 
            bridge[k]='X';
        }
    }
    free(bridge);
    return nkmers;
}

#endif
