#ifndef __ALPHABETS_HPP__
#define __ALPHABETS_HPP__

#include <string>
#include <vector>


namespace mtg {
namespace kmer {
namespace alphabets {

template <size_t n>
class log2 {
    static constexpr size_t _log2(size_t x) {
        if (x < 2) {
            return 0;
        } else {
            return _log2(x >> 1) + 1;
        }
    }

  public:
    static constexpr size_t value = _log2(n);
};


constexpr char kBOSSAlphabetProtein[] = "$ABCDEFGHIJKLMNOPQRSTUVWYZX";
constexpr uint8_t kBOSSSigmaProtein = sizeof(kBOSSAlphabetProtein) - 1;
constexpr uint8_t kBOSSBitsPerCharProtein = log2<kBOSSSigmaProtein - 1>::value + 1;
constexpr uint8_t kBOSSCharToProtein[128] = {
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
    16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26,
    26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
    16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26
};
const std::vector<uint8_t> kBOSSComplementMapProtein = {};
static_assert(kBOSSSigmaProtein <= 1llu << kBOSSBitsPerCharProtein);
static_assert(kBOSSSigmaProtein > 1llu << (kBOSSBitsPerCharProtein - 1));


//for case-specific DNA and RNA (U <-> T) data
constexpr char kBOSSAlphabetDNACaseSent[] = "$ACGTNacgt";
constexpr uint8_t kBOSSSigmaDNACaseSent = sizeof(kBOSSAlphabetDNACaseSent) - 1;
constexpr uint8_t kBOSSBitsPerCharDNACaseSent = log2<kBOSSSigmaDNACaseSent - 1>::value + 1;
constexpr uint8_t kBOSSCharToDNACaseSent[128] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 6, 5, 7,  5, 5, 5, 8,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  9, 9, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};
const std::vector<uint8_t> kBOSSComplementMapDNACaseSent = { 0, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
static_assert(kBOSSSigmaDNACaseSent <= 1llu << kBOSSBitsPerCharDNACaseSent);
static_assert(kBOSSSigmaDNACaseSent > 1llu << (kBOSSBitsPerCharDNACaseSent - 1));


//for DNA and RNA (U <-> T) alphabets
constexpr char kBOSSAlphabetDNA[] = "$ACGT";
constexpr uint8_t kBOSSSigmaDNA = sizeof(kBOSSAlphabetDNA) - 1;
constexpr uint8_t kBOSSBitsPerCharDNA = log2<kBOSSSigmaDNA - 1>::value + 1;
constexpr uint8_t kBOSSCharToDNA[128] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};
const std::vector<uint8_t> kBOSSComplementMapDNA = { 0, 4, 3, 2, 1, 5 };
static_assert(kBOSSSigmaDNA <= 1llu << kBOSSBitsPerCharDNA);
static_assert(kBOSSSigmaDNA > 1llu << (kBOSSBitsPerCharDNA - 1));

constexpr char kBOSSAlphabetDNA5[] = "$ACGTN";
constexpr uint8_t kBOSSSigmaDNA5 = sizeof(kBOSSAlphabetDNA5) - 1;
constexpr uint8_t kBOSSBitsPerCharDNA5 = log2<kBOSSSigmaDNA5 - 1>::value + 1;
static_assert(kBOSSSigmaDNA5 <= 1llu << kBOSSBitsPerCharDNA5);
static_assert(kBOSSSigmaDNA5 > 1llu << (kBOSSBitsPerCharDNA5 - 1));



constexpr char kAlphabetProtein[] = "ABCDEFGHIJKLMNOPQRSTUVWYZX";
constexpr uint8_t kSigmaProtein = sizeof(kAlphabetProtein) - 1;
constexpr uint8_t kBitsPerCharProtein = log2<kSigmaProtein - 1>::value + 1;
constexpr uint8_t kCharToProtein[128] = {
    25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,
    25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,
    25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,
    25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,  25, 25, 25, 25,
    25,  0,  1,  2,   3,  4,  5,  6,   7,  8,  9, 10,  11, 12, 13, 14,
    15, 16, 17, 18,  19, 20, 21, 22,  25, 23, 24, 25,  25, 25, 25, 25,
    25,  0,  1,  2,   3,  4,  5,  6,   7,  8,  9, 10,  11, 12, 13, 14,
    15, 16, 17, 18,  19, 20, 21, 22,  25, 23, 24, 25,  25, 25, 25, 25
};
const std::vector<uint8_t> kComplementMapProtein = {};
static_assert(kSigmaProtein <= 1llu << kBitsPerCharProtein);
static_assert(kSigmaProtein > 1llu << (kBitsPerCharProtein - 1));


//for case-specific DNA and RNA (U <-> T) data
constexpr char kAlphabetDNACaseSent[] = "ACGTNacgt";
constexpr uint8_t kSigmaDNACaseSent = sizeof(kAlphabetDNACaseSent) - 1;
constexpr uint8_t kBitsPerCharDNACaseSent = log2<kSigmaDNACaseSent - 1>::value + 1;
constexpr uint8_t kCharToDNACaseSent[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 5, 4, 6,  4, 4, 4, 7,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  8, 8, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
const std::vector<uint8_t> kComplementMapDNACaseSent = { 8, 7, 6, 5, 4, 3, 2, 1, 0 };
static_assert(kSigmaDNACaseSent <= 1llu << kBitsPerCharDNACaseSent);
static_assert(kSigmaDNACaseSent > 1llu << (kBitsPerCharDNACaseSent - 1));


//for DNA and RNA (U <-> T) alphabets
constexpr char kAlphabetDNA[] = "ACGT";
constexpr uint8_t kSigmaDNA = sizeof(kAlphabetDNA) - 1;
constexpr uint8_t kBitsPerCharDNA = log2<kSigmaDNA - 1>::value + 1;
constexpr uint8_t kCharToDNA[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
const std::vector<uint8_t> kComplementMapDNA = { 3, 2, 1, 0, 4 };
static_assert(kSigmaDNA <= 1llu << kBitsPerCharDNA);
static_assert(kSigmaDNA > 1llu << (kBitsPerCharDNA - 1));

constexpr char kAlphabetDNA5[] = "ACGTN";
constexpr uint8_t kSigmaDNA5 = sizeof(kAlphabetDNA5) - 1;
constexpr uint8_t kBitsPerCharDNA5 = log2<kSigmaDNA5 - 1>::value + 1;
static_assert(kSigmaDNA5 <= 1llu << kBitsPerCharDNA5);
static_assert(kSigmaDNA5 > 1llu << (kBitsPerCharDNA5 - 1));

} // namespace alphabets
} // namespace kmer
} // namespace mtg

#endif // __ALPHABETS_HPP__
