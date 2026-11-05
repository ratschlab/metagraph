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

// The rows a discovery decodes at once
constexpr size_t kResolveBatchRows = 4096;

struct ResolveOptions {
    // exactly one of |labels| (explicit) or |discover| must be set
    std::vector<std::string> labels;
    bool discover = false;
    size_t discover_max_labels = 1000;
    LabelKind discover_kind = LabelKind::COLUMN;
    // the shortest support run (in k-mers) that becomes a seed candidate
    uint64_t min_block_kmers = 1;
};

struct LabelProfile {
    LabelRef label;
    uint64_t kmers_supported = 0;
    std::vector<KmerInterval> runs;          // maximal support runs
};

struct SeedCandidate {
    KmerInterval kmers;
    std::vector<LabelId> labels;             // indices into SupportProfile::labels, sorted
};

// A discovery that found more labels than |discover_max_labels| keeps those supporting the
// most k-mers (ties: smaller column, then sequence id)
struct LabelTruncation {
    size_t kept = 0;
    size_t total = 0;
    uint64_t min_kept_kmers = 0;
    uint64_t max_dropped_kmers = 0;
    size_t dropped_full_length = 0;          // dropped labels supporting every in-graph k-mer
};

struct SupportProfile {
    size_t k = 0;
    uint64_t num_kmers = 0;
    Regime regime = Regime::BASIC;
    std::vector<KmerInterval> graph_runs;    // runs of k-mers present in the graph
    std::vector<LabelProfile> labels;
    std::optional<LabelTruncation> labels_truncated;
    std::vector<SeedCandidate> candidates;   // ordered: longer, more labels, smaller begin

    // base interval of a k-mer interval
    std::pair<uint64_t, uint64_t> bp_interval(const KmerInterval &iv) const {
        return { iv.begin, iv.end + k - 1 };
    }
};

// Resolve where |query| is supported by which labels. No traversal.
// Throws std::invalid_argument for unknown labels or invalid options.
SupportProfile resolve_support(const LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options);

struct SelectionPolicy {
    enum Policy { LONGEST_FIRST, MAX_SUPPORT };

    Policy policy = MAX_SUPPORT;
    size_t max_seeds = 10;
    uint64_t min_block_bp = 0;               // 0: at least k
    size_t max_labels_per_seed = 1000;
    bool merge_overlapping = true;
};

// The labels supporting a selected seed beyond |max_labels_per_seed|: the kept ones are the
// first in the order of their FNV-1a hash, so a cut is reproducible and unbiased by name
struct LabelPopulation {
    size_t supporting_total = 0;             // labels supporting the whole interval
    size_t included = 0;
    std::vector<std::string> dropped;        // kept only if <= 1000
    size_t dropped_count = 0;
    std::string dropped_digest;
};

struct FrozenSeed {
    std::string seed_id;
    std::string sequence;
    KmerInterval kmers;
    std::vector<std::string> labels;         // bytewise sorted
    LabelPopulation population;
    std::vector<size_t> overlaps_with;       // indices of other seeds sharing k-mers
};

struct SeedSelection {
    SelectionPolicy policy;
    size_t num_candidates = 0;
    size_t num_eligible = 0;
    std::vector<FrozenSeed> seeds;
};

SeedSelection select_seeds(const SupportProfile &profile,
                           std::string_view query,
                           const SelectionPolicy &policy,
                           bool canonical_orientation);

// seed_id = hex FNV-1a-64 over the case-mapped (and in canonical graphs canonical)
// sequence and the length-prefixed sorted labels
std::string make_seed_id(std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels);

uint64_t fnv1a64(std::string_view data, uint64_t hash = 0xcbf29ce484222325ULL);
// lower-case 16-digit hex of a 64-bit hash (seed ids, dropped-label digests)
std::string hex64(uint64_t x);

// Run-length encoding of a presence mask ("x<n>o<n>..."), compatible with the
// query output of `with_signature`.
std::string encode_runs(const std::vector<KmerInterval> &runs, uint64_t num_kmers);

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
