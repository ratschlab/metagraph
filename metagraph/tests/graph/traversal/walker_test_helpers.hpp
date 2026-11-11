#ifndef __TESTS_GRAPH_TRAVERSAL_WALKER_TEST_HELPERS_HPP__
#define __TESTS_GRAPH_TRAVERSAL_WALKER_TEST_HELPERS_HPP__

#include <algorithm>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "graph/traversal/walker.hpp"
#include "tests/graph/traversal/traversal_test_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace test {

constexpr size_t kK = 11;
constexpr size_t kLeft = static_cast<size_t>(graph::traversal::Arm::LEFT);
constexpr size_t kRight = static_cast<size_t>(graph::traversal::Arm::RIGHT);

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

// one sequence record (FASTA header and sequence) of an annotation column
struct Record {
    std::string column;
    std::string header;
    std::string sequence;
};

// A coordinate annotation of |records| with their headers (CoordToHeader), as built from
// FASTA files: a column's records, listed together, get consecutive coordinates
template <class Graph, class Annotation>
std::unique_ptr<graph::AnnotatedDBG>
build_records(size_t k, const std::vector<Record> &records,
              graph::DeBruijnGraph::Mode mode = graph::DeBruijnGraph::BASIC) {
    std::vector<std::string> sequences, columns;
    std::vector<uint64_t> starts;
    std::map<std::string, uint64_t> num_coords;
    for (const Record &r : records) {
        sequences.push_back(r.sequence);
        columns.push_back(r.column);
        starts.push_back(num_coords[r.column]);
        num_coords[r.column] += r.sequence.size() - k + 1;
    }
    auto anno = build_anno_graph<Graph, Annotation>(k, sequences, columns, mode, true, starts);
    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(num_coords.size());
    std::vector<std::vector<uint64_t>> num_kmers(num_coords.size());
    for (const Record &r : records) {
        const size_t column = encoder.encode(r.column);
        headers[column].push_back(r.header);
        num_kmers[column].push_back(r.sequence.size() - k + 1);
    }
    return with_coord_to_header(std::move(anno), std::move(headers), std::move(num_kmers));
}

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
