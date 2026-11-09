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
    // Not a request field: an explicit label the annotation does not hold is dropped and
    // listed in SupportProfile::unknown_labels (a multi-graph server, whose pairs need not hold
    // every label a request names: SPEC-multigraph-envelope.md section 2) instead of refusing
    // the request (std::invalid_argument "Label not found in the annotation"; a single-graph
    // server, the CLI). With every label absent the profile has no label
    bool drop_labels_not_in_graph = false;
    bool discover = false;
    // at least 1: a discovery keeps the top labels by k-mers supported, ties by name (bytes)
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
    // (DESIGN-traverse-graphlet.md §13): a prefix, never a sample of the whole query. Under it
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

    // The bounds of a discovery (SPEC-labeled-traversal-core.md §4.6; the server's
    // --resolve-max-discover-labels and --resolve-max-discover-memory-mb, the CLI none). Not
    // request fields.
    // |discover_max_tracked|: the most labels the discovery tracks (0: no cap). A label is
    // tracked from the first k-mer it is met at; the labels first met at a k-mer are admitted
    // all or none: the first k-mer whose new labels would make more than this many is the cut
    // k-mer, and no label first met there or later is tracked (SupportProfile::labels_capped)
    size_t discover_max_tracked = 0;
    // |discover_max_memory|: the bytes the discovery's pass may hold (0: no bound), charged
    // before they are allocated — the rows by the budget-aware reads, the rows kept for
    // repeats, the label index, the accumulators (kDiscover*Bytes) and the pass's arrays. The
    // pass decodes a batch in halves, gives up the rows kept and reads the next row alone
    // before it stops: only a k-mer whose row read alone, or whose accumulation, does not fit
    // beside the accumulators of the k-mers before it ends the work there (ResolveStop::MEMORY:
    // the resolve of the query's prefix before that k-mer)
    uint64_t discover_max_memory = 0;
};

// What a discovery charges against ResolveOptions::discover_max_memory for its accumulators
// (§4.6; the rows are charged by the decoder): per label tracked its accumulator and its id
// (held in deques: no relocation), and of a header label its entry in the label index (a hash
// table at most 90% full, which doubles: its entry, the empty buckets and the table before
// the growth while it is moved); per run and per trace break twice the element (a vector at
// most twice its size); per chain coordinate beyond those held inline its heap (charged as
// the vector's capacity grows). A model of the heap, never a measurement
constexpr uint64_t kDiscoverIndexEntryBytes = 112;
constexpr uint64_t kDiscoverRunBytes = 32;
constexpr uint64_t kDiscoverBreakBytes = 16;
// Under support "trace", the most a coordinate of a row adds to the chains' coordinates (8 B, a
// vector at most twice its size), which the check before a k-mer's accumulation counts with the
// growth of the scratch the coordinates are chained and mapped in (twice the row's largest
// column, 8 B a coordinate, each buffer). A presence discovery keeps no coordinate: only the
// sort scratch of an unsorted column (8 B a coordinate, charged before it grows)
constexpr uint64_t kDiscoverChainBytes = 16;

// The start of the refusal of an explicit seed naming a label the profile does not hold
// (select_seeds and the first-fault order of an explicit selection); the route adds to it what
// the answer-labels clamp left out (SPEC-labeled-traversal-core.md §4.6)
constexpr char kSeedLabelNotProfiled[] = "Explicit seed label was not profiled: '";

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

// A discovery's cut at ResolveOptions::discover_max_labels: the labels ranked by k-mers
// supported (more first), ties by name (bytes), the first |kept| of them kept. The k-mers
// supported are LabelProfile::kmers_supported's: under Support::TRACE the traced k-mers,
// else the k-mers the label is present at; every count below is on that scale
struct LabelTruncation {
    size_t kept = 0;
    size_t total = 0;
    uint64_t min_kept_kmers = 0;
    uint64_t max_dropped_kmers = 0;
    size_t dropped_full_length = 0;          // dropped labels supporting every in-graph k-mer
    // The labels with the last kept label's support (min_kept_kmers), kept or dropped, and how
    // many of them were kept: the tie group the cut went through (tied_kept < tied_at_boundary,
    // the kept ones the first by name) or ended with (tied_kept == tied_at_boundary)
    size_t tied_at_boundary = 0;
    size_t tied_kept = 0;
};

// A discovery's cut by ResolveOptions::discover_max_tracked (SPEC-labeled-traversal-core.md
// §4.6): the labels tracked are exactly those first met at a k-mer before |cut_kmer|, each
// counted over the whole query (its kmers_supported, runs and trace breaks are the uncut
// discovery's); a label first met at |cut_kmer| or later was not tracked, and supports at most
// |untracked_max_kmers| k-mers (the in-graph k-mers from |cut_kmer| on, of the prefix under a
// stop). LabelTruncation describes the tracked labels
struct LabelCap {
    size_t max_tracked = 0;
    size_t tracked = 0;
    uint64_t cut_kmer = 0;
    size_t new_at_cut = 0;                   // the labels first met at |cut_kmer|
    // the labels present at |cut_kmer| (new and tracked): above |max_tracked|, a discovery of
    // the sequence from |cut_kmer| meets more than the cap at its first k-mer and tracks nothing
    size_t labels_at_cut = 0;
    uint64_t untracked_max_kmers = 0;
    // the kept labels are the uncut discovery's top discover_max_labels: that many kept, the
    // last one supporting more than |untracked_max_kmers|
    bool ranking_exact = false;
    size_t max_labels = 0;                   // ResolveOptions::discover_max_labels
    size_t kept = 0;
    uint64_t min_kept_kmers = 0;             // of the last kept label (0: none kept)
    uint64_t met_at_least() const { return tracked + new_at_cut; }
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
    // TIME: the deadline (ResolveOptions::time_up); MEMORY: a discovery's memory bound
    // (ResolveOptions::discover_max_memory, |max_memory| bytes; always in phase ROWS)
    enum Reason { TIME, MEMORY };
    Reason reason = TIME;
    uint64_t max_memory = 0;
    // MEMORY at the query's first k-mer in the graph: the sequence from resolved_kmers has the
    // same k-mers in the graph, the same setup and the same first row, so a resolve of it stops
    // there again under the same bound (SPEC-labeled-traversal-core.md §4.6)
    bool cannot_progress = false;
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
    // the explicit labels the annotation does not hold, each once, in the request's order
    // (ResolveOptions::drop_labels_not_in_graph; empty otherwise)
    std::vector<std::string> unknown_labels;
    std::optional<LabelTruncation> labels_truncated;
    // set when the discovery's label cap bound (ResolveOptions::discover_max_tracked)
    std::optional<LabelCap> labels_capped;
    std::vector<SeedCandidate> candidates;   // ordered: longer, more labels, smaller begin
    // set when the deadline or a discovery's memory bound ended the work: every field above is
    // then exactly the resolve of the query's first num_kmers k-mers (a run ending there may
    // continue past it)
    std::optional<ResolveStop> stop;
    // The modelled memory of the work at its largest (bytes; work.memory_bytes, SPEC §4.4 and
    // §7.2.2): what the request's structures held at once — the query's k-mers and keys, the
    // decoded rows of a batch with the rows kept for later occurrences, a discovery's or the
    // explicit labels' accumulators, the explicit labels' primed hits and their fetch, and
    // the profile — each priced from its sizes and element sizes (profile_bytes,
    // selection_bytes), never from the allocator: it depends only on the index and the request
    // (and under a deadline on where the work stopped)
    uint64_t memory_peak = 0;

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

// What resolve_support refuses before it reads a row, in its order and words: the options
// (labels or discovery, a discovery keeping a label, trace on an index with coordinates in
// the BASIC regime), a query shorter than k, an explicit label the index does not hold (unless
// options.drop_labels_not_in_graph), a header discovery without a record mapping. Throws
// std::invalid_argument; reads nothing of the matrix (a single-graph /resolve is admitted
// after it)
void check_resolve_request(const LabelOracle &oracle,
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

// What a profile and a selection hold, as SupportProfile::memory_peak prices them: their
// fixed parts, every label's name and runs, the candidates' label lists; a seed's id, sequence
// and label names
uint64_t profile_bytes(const SupportProfile &profile);
uint64_t selection_bytes(const SeedSelection &selection);

SeedSelection select_seeds(const SupportProfile &profile,
                           std::string_view query,
                           const SelectionPolicy &policy,
                           bool canonical_orientation);

// What select_seeds refuses of an explicit selection that |query|'s length, k, |options|' labels
// and |policy| decide alone: a seed's kmer_interval empty or out of range, a seed's label that
// is none of the explicit labels (nothing with a discovery's labels, known after the work, or
// with options.drop_labels_not_in_graph). Throws std::invalid_argument in select_seeds' words,
// and only when a refusal is certain: the first fault in select_seeds' order is the one stated,
// which for a request with two faults may be a seed not fully in the graph (found by mapping
// the query's k-mers; no row is read). A seed not fully in the graph is otherwise left to
// select_seeds. Call after check_resolve_request (|query| is not shorter than k)
void check_explicit_seeds(const LabelOracle &oracle,
                          std::string_view query,
                          const ResolveOptions &options,
                          const SelectionPolicy &policy);

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
