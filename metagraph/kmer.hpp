#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <immintrin.h>
#include <sstream>

#include <htslib/kseq.h>

typedef uint64_t TAlphabet;

//NOTE: assumes that kBitsPerChar <= 64
const int kBitsPerChar = 3;
const int kMax = 1llu << kBitsPerChar;

class KMer {
    friend std::ostream& operator<<(std::ostream &os, const KMer &kmer);

    template <typename T>
    friend KMer& operator<<(KMer &kmer, const T &shift);
  public:
    template <typename T>
    KMer(const T &arr, size_t k) {
        memset(s_, 0, sizeof(s_));
        assert(k * kBitsPerChar < 256 && k >= 2
                && "String must be between lengths 2 and 256 / kBitsPerChar");
   
        for (int i = k - 2; i >= 0; --i) {
            assert(arr[i] + 1 < kMax && "Alphabet size too big for the given number of bits");
            *this <<= kBitsPerChar;
            s_[0] |= arr[i] + 1;
        }
        if (k >= 1) {
            *this <<= kBitsPerChar;
            s_[0] |= arr[k - 1] + 1;
        }
    }

    template <class Map, class String>
    KMer(const String &seq, Map &&to_alphabet);

    KMer(KMer &&other) {
        memcpy(s_, other.s_, sizeof(other.s_));
    }
    KMer(const KMer &other) {
        memcpy(s_, other.s_, sizeof(other.s_));
    }

    explicit KMer(const size_t &a, const size_t &b, const size_t &c, const size_t &d) {
        //seq_ = _mm256_set_epi64x(d,c,b,a);
        s_[0] = a;
        s_[1] = b;
        s_[2] = c;
        s_[3] = d;
    }

    explicit KMer(const size_t &a) {
        s_[0] = a;
        s_[1] = 0;
        s_[2] = 0;
        s_[3] = 0;
    }

    const size_t& nth(const size_t &i) const {
        return s_[i];
    }

    size_t& nth(const size_t &i) {
        return s_[i];
    }


    KMer& operator=(KMer &&other) {
        if (&(s_[0]) != &(other.s_[0])) {
            //TODO: vectorize this
            memcpy(s_, other.s_, sizeof(s_));
        }
        return *this;
    }
    KMer& operator=(const KMer &other) {
        if (&(s_[0]) != &(other.s_[0])) {
            //TODO: vectorize this
            memcpy(s_, other.s_, sizeof(s_));
        }
        return *this;
    }

    KMer operator~() {
        KMer that(*this);
        that.flip();
        return that;
    }

    KMer& flip() {
        for (uint8_t i = 0; i < 4; ++i) {
            s_[i] = ~s_[i];
        }
        return *this;
    }

    bool operator<(const KMer &other) const {
        int i;
        for (i = 3; i >= 0 && s_[i] == other.s_[i]; --i);
        return (i >= 0) && (s_[i] < other.s_[i]);
    }
    bool operator!=(const KMer &other) const {
        //TODO: vectorize this
        for (int i = 3; i >= 0; --i) {
            if (s_[i] != other.s_[i])
                return 1;
        }
        return 0;
    }
    bool operator==(const KMer &other) const {
        return !(*this != other);
    }

    TAlphabet operator[](size_t i) const { return get(i) - 1; }

    static bool compare_kmer_suffix(const KMer &k1,
                                    const KMer &k2, size_t minus = 0);

    std::string to_string(const std::string &alphabet) const;

    std::string to_ints() const {
        std::ostringstream os;
        for (uint8_t i = 0; i < 4; ++i) {
            os << s_[i] << " ";
        }
        return os.str();
    }

    //these functions piece together blocks when indexing between them
    KMer& operator<<=(const int shift) {
        if (shift >= 256) {
            memset(s_, 0, sizeof(s_));
        } else if (shift % 64 == 0) {
            int i = 3;
            for (; i >= shift / 64; --i) {
                s_[i] = s_[i - shift / 64];
            }
            for (; i >= 0; --i) {
                s_[i] = 0;
            }
        } else {
            int i = 3;
            for (; i > shift / 64; --i) {
                s_[i] = (s_[i - shift / 64] & ((1llu << (64 - (shift % 64))) - 1));
                s_[i] <<= (shift % 64);
                s_[i] |= s_[i - shift / 64 - 1] >> (64 - (shift % 64));
            }
            s_[i] = (s_[i - shift / 64] & ((1llu << (64 - (shift % 64))) - 1)) << (shift % 64);
            for (--i; i >= 0; --i) {
                s_[i] = 0;
            }
        }
        return *this;
    }


    KMer& operator>>=(const int shift) {
        if (shift >= 256) {
            memset(s_, 0, sizeof(s_));
        } else if (shift % 64 == 0) {
            int i = 0;
            for (; i + shift / 64 < 4; ++i) {
                s_[i] = s_[i + shift / 64];
            }
            for (; i < 4; ++i) {
                s_[i] = 0;
            }
        } else {
            int i = 0;
            for (; i + shift / 64 < 3; ++i) {
                s_[i] = s_[i + shift / 64] >> (shift % 64);
                s_[i] |= (s_[i + shift / 64 + 1] & ((1llu << (shift % 64)) - 1)) << (64 - (shift % 64));
            }
            s_[i] = s_[i + shift / 64] >> (shift % 64);
            for (++i; i < 4; ++i) {
                s_[i] = 0;
            }
        }
        return *this;
    }

    KMer& operator&=(const KMer &other) {
        for (uint8_t i = 0; i < 4; ++i) {
            s_[i] &= other.s_[i];
        }
        return *this;
    }

    KMer& operator|=(const KMer &other) {
        for (uint8_t i = 0; i < 4; ++i) {
            s_[i] |= other.s_[i];
        }
        return *this;
    }

    //When working with int types
    KMer& operator&=(const size_t &other) {
        s_[0] &= other;
        for (uint8_t i = 1; i < 4; ++i) {
            s_[i] = 0;
        }
        return *this;
    }
    template<typename T>
    KMer& operator|=(const T &other) {
        s_[0] |= other;
        return *this;
    }
    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     */

    template<typename T>
    KMer& update(const size_t &k, const T &edge) {
        KMer last(this->get(0));
        last <<= kBitsPerChar * k;
        *this >>= kBitsPerChar;
        *this |= last;
        this->nth(0) &= kErasingFirstCodeMask;
        *this |= edge + 1;
        return *this;
    }

  private:
    uint64_t s_[4]; // kmer sequence
    //__m256i seq_;

    TAlphabet get(size_t i) const;

    static const long kErasingFirstCodeMask = ~((1l << kBitsPerChar) - 1);

    };


template <class Map, class String>
KMer::KMer(const String &seq, Map &&to_alphabet) {
    std::vector<uint8_t> arr(seq.size());
    std::transform(seq.begin(), seq.end(), arr.begin(), to_alphabet);
    *this = KMer(arr.data(), arr.size());
}

#endif // __KMER_HPP__
