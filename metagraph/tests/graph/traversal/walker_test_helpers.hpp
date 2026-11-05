#ifndef __TESTS_GRAPH_TRAVERSAL_WALKER_TEST_HELPERS_HPP__
#define __TESTS_GRAPH_TRAVERSAL_WALKER_TEST_HELPERS_HPP__

#include <algorithm>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "graph/traversal/walker.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace test {

constexpr size_t kK = 11;
constexpr size_t kLeft = static_cast<size_t>(graph::traversal::Arm::LEFT);
constexpr size_t kRight = static_cast<size_t>(graph::traversal::Arm::RIGHT);

inline std::string rc(std::string s) {
    ::reverse_complement(s);
    return s;
}

inline std::string random_seq(size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(len, 'A');
    for (char &c : s) {
        c = "ACGT"[gen() % 4];
    }
    return s;
}

// No repeated (k-1)-mer or k-mer in either orientation and no reverse-complement
// palindromes of length k-1 or k+1: every k-mer is one node with one successor in every
// graph mode, and no hairpin
inline bool is_clean(const std::string &s) {
    for (size_t len : { kK - 1, kK }) {
        std::set<std::string> seen;
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string km = s.substr(i, len);
            if (seen.count(km) || seen.count(rc(km)))
                return false;
            seen.insert(km);
        }
    }
    for (size_t len : { kK - 1, kK + 1 }) {
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string x = s.substr(i, len);
            if (rc(x) == x)
                return false;
        }
    }
    return true;
}

inline std::string clean_block(size_t len, uint32_t seed) {
    for (uint32_t s = seed; ; ++s) {
        std::string b = random_seq(len, s);
        if (is_clean(b))
            return b;
    }
}

// blocks of the given lengths cut from one clean sequence: clean together, in any order
// that keeps their junctions apart
inline std::vector<std::string> clean_blocks(const std::vector<size_t> &lengths, uint32_t seed) {
    size_t total = 0;
    for (size_t l : lengths) {
        total += l;
    }
    std::string master = clean_block(total, seed);
    std::vector<std::string> blocks;
    size_t pos = 0;
    for (size_t l : lengths) {
        blocks.push_back(master.substr(pos, l));
        pos += l;
    }
    return blocks;
}

// X, P, Q with P and Q starting with different bases: X + P and X + Q fork after X
inline std::vector<std::string> fork_blocks(uint32_t seed,
                                            std::vector<size_t> lengths = { 30, 40, 40 }) {
    for (uint32_t s = seed; ; ++s) {
        auto b = clean_blocks(lengths, s);
        if (b[1][0] != b[2][0])
            return b;
    }
}

// A on W·R·X, B on R·X·Y
struct LabelEndFixture {
    std::string W, R, X, Y;
    LabelEndFixture() {
        auto b = clean_blocks({ 20, 60, 25, 30 }, 2);
        W = b[0]; R = b[1]; X = b[2]; Y = b[3];
    }
    std::vector<std::string> seqs() const { return { W + R + X, R + X + Y }; }
    std::vector<std::string> labels() const { return { "A", "B" }; }
};

// ... and C on the first 40 bases of R
struct DerivedFixture : public LabelEndFixture {
    std::string prefix() const { return R.substr(0, 40); }
    std::vector<std::string> seqs() const { return { W + R + X, R + X + Y, prefix() }; }
    std::vector<std::string> labels() const { return { "A", "B", "C" }; }
};

// the bases of a path from the seed boundary outward (natural orientation on the left arm)
inline std::string spell_path(const graph::traversal::ArmResult &arm,
                              const graph::traversal::PathResult &path) {
    std::vector<const std::string*> pieces;
    for (size_t s = path.leaf; ; s = arm.segments[s].parents[0]) {
        pieces.push_back(&arm.segments[s].sequence);
        if (arm.segments[s].parents.empty())
            break;
    }
    if (arm.arm == graph::traversal::Arm::RIGHT)
        std::reverse(pieces.begin(), pieces.end());
    std::string out;
    for (const std::string *p : pieces) {
        out += *p;
    }
    return out;
}

} // namespace test
} // namespace mtg

#endif // __TESTS_GRAPH_TRAVERSAL_WALKER_TEST_HELPERS_HPP__
