#ifndef __TRAVERSAL_RESOLVE_HPP__
#define __TRAVERSAL_RESOLVE_HPP__

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace mtg::graph { class AnnotatedDBG; }
namespace mtg::graph::traversal {

// Extracted from gr/labeled-traversal, 07dfa4875fcf. Intervals are zero-based,
// half-open k-mer-start intervals; membership does not establish a record trace.
using LabelId = uint32_t;
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


struct LabelRef { uint64_t column = 0; std::string name; };
struct LabelProfile {
    LabelRef label;
    uint64_t kmers_supported = 0;
    std::vector<KmerInterval> runs;
};
struct SeedCandidate { KmerInterval kmers; std::vector<LabelId> labels; };
struct LabelTruncation {
    size_t kept = 0, total = 0;
    uint64_t min_kept_kmers = 0, max_dropped_kmers = 0;
    size_t dropped_full_length = 0;
};
struct SupportProfile {
    size_t k = 0;
    uint64_t num_kmers = 0;
    std::vector<KmerInterval> graph_runs;
    std::vector<LabelProfile> labels;
    std::optional<LabelTruncation> labels_truncated;
    std::vector<SeedCandidate> candidates;
    std::pair<uint64_t, uint64_t> bp_interval(const KmerInterval &iv) const {
        return { iv.begin, iv.end + k - 1 };
    }
};
struct ResolveOptions {
    std::vector<std::string> labels;
    bool discover = false;
    size_t discover_max_labels = 1000;
    uint64_t min_block_kmers = 1;
};

// Profile annotation columns, independently of an optional .seqs mapping.
// Discovery ranks by matched k-mers, ties by column id, and reports truncation.
// Limits on returned labels are not limits on annotation decoding or memory.
SupportProfile resolve_support(const AnnotatedDBG &, std::string_view,
                               const ResolveOptions &);

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


std::string make_seed_id(const std::string &release_id, std::string_view sequence,
                         bool canonical_orientation, std::vector<std::string> labels);
std::string encode_runs(const std::vector<KmerInterval> &, uint64_t num_kmers);

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

} // namespace mtg::graph::traversal
#endif
