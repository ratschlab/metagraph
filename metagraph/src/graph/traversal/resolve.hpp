#ifndef __TRAVERSAL_RESOLVE_HPP__
#define __TRAVERSAL_RESOLVE_HPP__

#include <optional>
#include <string>
#include <string_view>
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

    Support support = Support::KMER;
    uint64_t min_block_kmers = 1;
};

struct LabelProfile {
    LabelRef label;
    uint64_t kmers_supported = 0;
    std::vector<KmerInterval> runs;          // maximal support runs
    std::vector<uint64_t> trace_breaks;      // k-mer starts where a trace run ended early
};

struct SeedCandidate {
    KmerInterval kmers;
    std::vector<LabelId> labels;             // indices into SupportProfile::labels, sorted
};

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
    Support support = Support::KMER;
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
// Throws std::invalid_argument for unknown labels, invalid options, or unsupported
// support kinds (trace needs coordinates and the BASIC regime).
SupportProfile resolve_support(LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options);


struct ExplicitSeed {
    KmerInterval kmers;
    std::vector<std::string> labels;
};

struct SelectionPolicy {
    enum Policy { LONGEST_FIRST, MAX_SUPPORT, EXPLICIT };
    enum LabelOrder { HASH, COLUMN_ID, KMERS_SUPPORTED };

    Policy policy = MAX_SUPPORT;
    size_t max_seeds = 10;
    uint64_t min_block_bp = 0;               // 0: at least k
    size_t max_labels_per_seed = 1000;
    LabelOrder label_order = HASH;
    uint64_t sample_seed = 0;
    bool merge_overlapping = true;
    std::vector<ExplicitSeed> explicit_seeds;
    std::string release_id;                  // folded into seed ids
};

struct LabelPopulation {
    size_t supporting_total = 0;             // labels supporting the whole interval
    size_t included = 0;
    SelectionPolicy::LabelOrder order = SelectionPolicy::HASH;
    uint64_t sample_seed = 0;
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
    std::vector<std::string> labels_not_covering;  // EXPLICIT only
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

// seed_id = hex FNV-1a-64 over release, the case-mapped canonical sequence and the
// length-prefixed sorted labels
std::string make_seed_id(const std::string &release_id,
                         std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels);

uint64_t fnv1a64(std::string_view data, uint64_t hash = 0xcbf29ce484222325ULL);
static constexpr uint64_t kFnvOffsetBasis = 0xcbf29ce484222325ULL;
// lower-case 16-digit hex of a 64-bit hash (seed ids, dropped-label digests)
std::string hex64(uint64_t x);

// Run-length encoding of a presence mask ("x<n>o<n>..."), compatible with the
// query output of `with_signature` (moved here from cli/query.cpp).
std::string encode_runs(const std::vector<KmerInterval> &runs, uint64_t num_kmers);

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_RESOLVE_HPP__
