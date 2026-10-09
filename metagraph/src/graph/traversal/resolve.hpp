#ifndef __TRAVERSAL_RESOLVE_HPP__
#define __TRAVERSAL_RESOLVE_HPP__

#include <functional>
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

// A /resolve decodes the present k-mers' rows in batches and holds one batch at a time, so
// that the rows it holds stay bounded. The first batch has kResolveFirstBatchRows rows; each
// next one is sized from the widest row of the one before to about kResolveBatchBytes
// (row_copy_bytes), at most twice the previous batch's rows and at most kResolveBatchRows rows
// — so a batch is larger than kResolveBatchBytes only where the rows widen within the query. A
// row is decoded once per batch; a discovery keeps the row of a k-mer that occurs again later
// in the query for that occurrence, at most kResolveKeptBytes of such rows (beyond that a
// repeated row is decoded again), and an explicit profile decodes each distinct row once (its
// hits are kept per key)
constexpr size_t kResolveFirstBatchRows = 64;
constexpr size_t kResolveBatchRows = 4096;
constexpr uint64_t kResolveBatchBytes = uint64_t(64) << 20;
constexpr uint64_t kResolveKeptBytes = uint64_t(256) << 20;

struct ResolveOptions {
    // exactly one of |labels| (explicit) or |discover| must be set
    std::vector<std::string> labels;
    bool discover = false;
    size_t discover_max_labels = 1000;
    LabelKind discover_kind = LabelKind::COLUMN;

    Support support = Support::KMER;
    uint64_t min_block_kmers = 1;

    // Polled between the phases of resolve_support() (before and after the discovery read,
    // before and after the support fetch): true abandons the request, throwing |abandon|'s
    // exception — the server's check that the client is still connected. Not a request field.
    std::function<bool()> stop;
    std::function<void()> abandon;
    // The work deadline of a request with bounds.time_budget_ms (unset: no deadline, and the
    // request runs exactly as without one). Read where |stop| is polled during the work —
    // between two row batches, every kResolveCheckKmers k-mers of the explicit labels' support
    // pass —, after |stop|: true ends the work there, and the profile is then exactly the
    // resolve of the query's first SupportProfile::stop->resolved_kmers k-mers
    // (DESIGN-traverse-graphlet.md §21): a prefix, never a sample of the whole query. Under it
    // the explicit labels' hits are fetched kResolveCheckKmers k-mers at a time (one fetch of
    // the whole query, the unbudgeted path, is a piece no clock read can end). Not a request
    // field.
    std::function<bool()> time_up;
    // Read every kResolveCheckLabels labels of the loops after the work (a discovery's ranking
    // and naming, the profiles, the candidates' grouping): throws to abandon a request whose
    // answer can no longer be built and written in time (the server: 503 deadline). Not a
    // request field.
    std::function<void()> finish_check;
    // the batches the rows are decoded in: at most |batch_rows| rows, sized to about
    // |batch_bytes|, and the bound of the rows a discovery keeps for repeated k-mers (tests
    // vary them; the profile does not depend on them). Not request fields
    size_t batch_rows = kResolveBatchRows;
    uint64_t batch_bytes = kResolveBatchBytes;
    uint64_t kept_bytes = kResolveKeptBytes;
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

// The k-mers between two checkpoints of the explicit labels' support pass: a hit list per
// k-mer is a few to a few thousand hits, so a gone client (and a deadline) is seen within
// milliseconds, and a check (a peek on its socket, a clock read) costs nothing beside that
constexpr uint64_t kResolveCheckKmers = 4096;

// The loops over the labels after a /resolve's work read ResolveOptions::finish_check once in
// this many labels: a label's step there is a few allocations (its name, its runs), so the
// answer's time is read within a few milliseconds of it
constexpr size_t kResolveCheckLabels = 4096;

// Where a deadline (ResolveOptions::time_up) ended a /resolve's work, and how far it got: the
// profile is that of the query's first |resolved_kmers| k-mers (SupportProfile::num_kmers),
// out of its |query_kmers|
struct ResolveStop {
    // ROWS: reading the present k-mers' rows (a discovery's pass, the explicit labels'
    // priming) or before it; SUPPORT: the explicit labels' hits, k-mer by k-mer
    enum Phase { ROWS, SUPPORT };
    Phase phase = ROWS;
    uint64_t resolved_kmers = 0;
    uint64_t query_kmers = 0;
    // the runs of the WHOLE query's k-mers present in the graph (SupportProfile::graph_runs
    // holds the prefix's): the k-mers are mapped before the deadline is first read, so their
    // presence is known whatever the stop — an explicit seed not fully in the graph is a 400
    // under a stop as without one
    std::vector<KmerInterval> query_graph_runs;
};

struct SupportProfile {
    size_t k = 0;
    // the k-mers profiled: the query's, or under a stop its prefix's (stop->resolved_kmers)
    uint64_t num_kmers = 0;
    Regime regime = Regime::BASIC;
    Support support = Support::KMER;
    std::vector<KmerInterval> graph_runs;    // runs of k-mers present in the graph
    std::vector<LabelProfile> labels;
    std::optional<LabelTruncation> labels_truncated;
    std::vector<SeedCandidate> candidates;   // ordered: longer, more labels, smaller begin
    // set when the deadline ended the work: every field above is then exactly the resolve of
    // the query's first num_kmers k-mers (a run ending there may continue past it)
    std::optional<ResolveStop> stop;

    // base interval of a k-mer interval
    std::pair<uint64_t, uint64_t> bp_interval(const KmerInterval &iv) const {
        return { iv.begin, iv.end + k - 1 };
    }
};

// Resolve where |query| is supported by which labels. No traversal.
// Throws std::invalid_argument for unknown labels, invalid options, or unsupported
// support kinds (trace needs coordinates and the BASIC regime) — whatever the deadline: the
// labels are resolved before it is first read.
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
// query output of `with_signature`.
std::string encode_runs(const std::vector<KmerInterval> &runs, uint64_t num_kmers);

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_RESOLVE_HPP__
