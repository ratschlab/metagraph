#ifndef __TRAVERSAL_RESOLVE_HPP__
#define __TRAVERSAL_RESOLVE_HPP__

#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "traversal_types.hpp"
#include "label_oracle.hpp"


namespace mtg {
namespace graph {
namespace traversal {

// half-open interval of k-mer starts [begin, end)
struct KmerInterval {
    uint64_t begin = 0;
    uint64_t end = 0;

    uint64_t size() const { return end - begin; }
    bool operator==(const KmerInterval &o) const { return begin == o.begin && end == o.end; }
    bool operator<(const KmerInterval &o) const {
        return std::tie(begin, end) < std::tie(o.begin, o.end);
    }
};

struct ResolveOptions {
    // exactly one of |labels| (explicit) or |discover| must be set
    std::vector<std::string> labels;
    bool discover = false;
    size_t discover_max_labels = 1000;
    LabelKind discover_kind = LabelKind::COLUMN;
};

struct LabelProfile {
    LabelRef label;
    uint64_t kmers_supported = 0;
    std::vector<KmerInterval> runs;          // maximal support runs
};

// A discovery that found more labels than |discover_max_labels| keeps those supporting the
// most k-mers (ties: smaller column, then sequence id)
struct LabelTruncation {
    size_t kept = 0;
    size_t total = 0;
};

struct SupportProfile {
    size_t k = 0;
    uint64_t num_kmers = 0;
    Regime regime = Regime::BASIC;
    std::vector<KmerInterval> graph_runs;    // runs of k-mers present in the graph
    std::vector<LabelProfile> labels;
    std::optional<LabelTruncation> labels_truncated;
};

// Resolve where |query| is supported by which labels. No traversal.
// Throws std::invalid_argument for unknown labels or invalid options.
SupportProfile resolve_support(const LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options);

// Seeds cover the labels: the first is the interval contained in the most label runs, each
// next one the interval contained in the most runs of labels that no earlier seed carries
struct SelectionPolicy {
    size_t max_seeds = 10;
    uint64_t min_block_bp = 0;               // 0: at least k
    size_t max_labels_per_seed = 1000;
};

// The labels supporting a selected seed beyond |max_labels_per_seed|: the kept ones are the
// first in the order of their FNV-1a hash, so a cut is reproducible and unbiased by name
struct LabelPopulation {
    size_t supporting_total = 0;             // labels supporting the whole interval
    std::vector<std::string> dropped;        // kept only if <= 1000
    size_t dropped_count = 0;
};

struct FrozenSeed {
    std::string seed_id;
    std::string sequence;
    KmerInterval kmers;
    std::vector<std::string> labels;         // bytewise sorted
    LabelPopulation population;
};

// Seeds are listed in the order they were picked; they may overlap
std::vector<FrozenSeed> select_seeds(const SupportProfile &profile,
                                     std::string_view query,
                                     const SelectionPolicy &policy,
                                     bool canonical_orientation);

// seed_id = hex FNV-1a-64 over the case-mapped (and in canonical graphs canonical)
// sequence and the length-prefixed sorted labels
std::string make_seed_id(std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels);

// |sequence| in the case the graph stores (upper case unless the build is case-sensitive)
std::string in_graph_case(std::string_view sequence);

// The maximal runs of the positions i in [0, n) where |present|(i) holds, ascending
template <class Present>
std::vector<KmerInterval> runs_of(uint64_t n, const Present &present) {
    std::vector<KmerInterval> runs;
    for (uint64_t i = 0; i < n; ++i) {
        if (!present(i))
            continue;
        if (runs.empty() || runs.back().end != i) {
            runs.push_back({ i, i + 1 });
        } else {
            runs.back().end = i + 1;
        }
    }
    return runs;
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_RESOLVE_HPP__
