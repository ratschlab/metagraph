#ifndef __ALPHABETS_HPP__
#define __ALPHABETS_HPP__

#include <string>
#include <vector>


namespace alphabets {

constexpr size_t log2(size_t n) {
    if (n < 2) {
        return 0;
    } else {
        return log2(n >> 2) + 1;
    }
}


constexpr char kBOSSAlphabetProtein[] = "$ABCDEFGHIJKLMNOPQRSTUVWYZX";
constexpr uint8_t kBOSSLogSigmaProtein = log2(sizeof(kBOSSAlphabetProtein) - 2) + 1;
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
const std::vector<uint8_t> kBOSSCanonicalMapProtein = {};


//for case-specific DNA and RNA (U <-> T) data
constexpr char kBOSSAlphabetDNACaseSent[] = "$ACGTNacgt";
constexpr uint8_t kBOSSLogSigmaDNACaseSent = log2(sizeof(kBOSSAlphabetDNACaseSent) - 2) + 1;
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
const std::vector<uint8_t> kBOSSCanonicalMapDNACaseSent = { 0, 9, 8, 7, 6, 5, 4, 3, 2, 1 };


//for DNA and RNA (U <-> T) alphabets
constexpr char kBOSSAlphabetDNA[] = "$ACGTN";
constexpr uint8_t kBOSSLogSigmaDNA = log2(sizeof(kBOSSAlphabetDNA) - 2) + 1;
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
const std::vector<uint8_t> kBOSSCanonicalMapDNA = { 0, 4, 3, 2, 1, 5 };


//for DNA and RNA (U <-> T) alphabets
constexpr char kBOSSAlphabetDNA4[] = "$ACGT";
constexpr uint8_t kBOSSLogSigmaDNA4 = log2(sizeof(kBOSSAlphabetDNA4) - 2) + 1;
constexpr uint8_t kBOSSCharToDNA4[128] = {
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 2,  1, 1, 1, 3,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  4, 4, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 2,  1, 1, 1, 3,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  4, 4, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1
};
const std::vector<uint8_t> kBOSSCanonicalMapDNA4 = { 0, 4, 3, 2, 1 };



constexpr char kAlphabetProtein[] = "ABCDEFGHIJKLMNOPQRSTUVWYZX";
constexpr uint8_t kLogSigmaProtein = log2(sizeof(kAlphabetProtein) - 2) + 1;
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
const std::vector<uint8_t> kCanonicalMapProtein = {};


//for case-specific DNA and RNA (U <-> T) data
constexpr char kAlphabetDNACaseSent[] = "ACGTNacgt";
constexpr uint8_t kLogSigmaDNACaseSent = log2(sizeof(kAlphabetDNACaseSent) - 2) + 1;
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
const std::vector<uint8_t> kCanonicalMapDNACaseSent = { 8, 7, 6, 5, 4, 3, 2, 1, 0 };


//for DNA and RNA (U <-> T) alphabets
constexpr char kAlphabetDNA[] = "ACGTN";
constexpr uint8_t kLogSigmaDNA = log2(sizeof(kAlphabetDNA) - 2) + 1;
constexpr uint8_t kCharToDNA[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
const std::vector<uint8_t> kCanonicalMapDNA = { 3, 2, 1, 0, 4 };


//for DNA and RNA (U <-> T) alphabets
constexpr char kAlphabetDNA4[] = "ACGT";
constexpr uint8_t kLogSigmaDNA4 = log2(sizeof(kAlphabetDNA4) - 2) + 1;
constexpr uint8_t kCharToDNA4[128] = {
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};
const std::vector<uint8_t> kCanonicalMapDNA4 = { 3, 2, 1, 0 };


} // namespace alphabets

#endif // __ALPHABETS_HPP__
