#include "resolve.hpp"

#include <algorithm>
#include <iterator>
#include <cctype>
#include <deque>
#include <map>
#include <numeric>
#include <optional>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "graph/annotated_dbg.hpp"
#include "graph/traversal/kept_rows.hpp"
#include "graph/traversal/seeded_hash.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_diff/row_diff_cache.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/logger.hpp"
#include "common/utils/string_utils.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;
using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;
using annot::matrix::DecodeBudget;
using annot::matrix::DecodeStatus;

uint64_t fnv1a64(std::string_view data, uint64_t hash) {
    for (unsigned char c : data) {
        hash ^= c;
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

std::string hex64(uint64_t x) {
    static const char *digits = "0123456789abcdef";
    std::string s(16, '0');
    for (int i = 15; i >= 0; --i) {
        s[i] = digits[x & 15];
        x >>= 4;
    }
    return s;
}

std::string make_seed_id(const std::string &release_id,
                         std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels) {
    // defined over the sequence after the build's case mapping (§4.3 / §6.1), so a
    // seed resubmitted in another case keeps its id
    std::string seq(sequence);
#if ! _DNA_CASE_SENSITIVE_GRAPH
    for (char &c : seq) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
#endif
    if (canonical_orientation) {
        std::string rc(seq);
        ::reverse_complement(rc);
        if (rc < seq)
            seq = rc;
    }
    std::sort(labels.begin(), labels.end());
    uint64_t h = fnv1a64(release_id);
    h = fnv1a64("\t", h);
    h = fnv1a64(seq, h);
    h = fnv1a64("\t", h);
    for (const auto &l : labels) {
        h = fnv1a64(std::to_string(l.size()) + ":" + l, h);
    }
    return hex64(h);
}

std::string encode_runs(const std::vector<KmerInterval> &runs, uint64_t num_kmers) {
    std::string out;
    uint64_t pos = 0;
    for (const auto &run : runs) {
        if (run.begin > pos)
            out += "o" + std::to_string(run.begin - pos);
        out += "x" + std::to_string(run.size());
        pos = run.end;
    }
    if (num_kmers > pos)
        out += "o" + std::to_string(num_kmers - pos);
    return out;
}

static std::vector<KmerInterval> runs_of(const std::vector<bool> &mask) {
    return runs_of(mask.size(), [&](uint64_t i) { return mask[i]; });
}


namespace {

/**
 * The support of one label, built k-mer by k-mer in ascending k-mer order (each k-mer at most
 * once): what a discovery accumulates while it reads the rows, and the trace of an explicit
 * label over its hits.
 *   add(i)                 the label is present at k-mer i: |kmers|, and the runs its maximal
 *                          runs (runs_of over the k-mers it is present at)
 *   add_trace(i, coords)   the label is present at k-mer i with |coords| (sorted; possibly
 *                          none): |kmers| as add, and the runs, |trace_breaks| and |traced| the
 *                          trace-consistent runs — a chain of coordinates increasing by one per
 *                          k-mer, a k-mer without coordinates (or without the label: a gap in
 *                          the calls) breaking it. |scratch|: a buffer the caller reuses (a
 *                          discovery calls this for every label of every row)
 *   take_runs()            the runs, once the last k-mer was added
 * A discovery calls these for every label of every row, so what a call reads and writes — the
 * counts, the run being extended and a chain of up to four coordinates — is kept together in
 * the accumulator, and its vectors are touched only when a run closes or a trace jumps
 */
struct LabelSupport {
    static constexpr uint64_t kNone = std::numeric_limits<uint64_t>::max();

    uint64_t kmers = 0;
    uint64_t traced = 0;
    uint64_t last = kNone;                   // the k-mer of the last add_trace
    KmerInterval open { kNone, kNone };      // the run being extended, not yet in |runs|
    SmallVector<Coord, 4> live;              // the chain's coordinates at k-mer |last|
    std::vector<KmerInterval> runs;          // the closed runs
    std::vector<uint64_t> trace_breaks;      // k-mer starts where a trace run ended early

    void add(uint64_t i) {
        kmers++;
        extend(i);
    }

    void add_trace(uint64_t i, const Coord *coords, size_t n, std::vector<Coord> &scratch) {
        kmers++;
        // the k-mers since the last call did not have the label: the chain ended there
        if (last == kNone || last + 1 != i)
            live.clear();
        last = i;
        if (!n) {
            live.clear();
            return;
        }
        if (live.size() == n) {
            // the common case on a conserved stretch: every coordinate continues the chain,
            // which moves on by one in place (what the merge below gives, without its copies)
            size_t j = 0;
            while (j < n && coords[j] == live[j] + 1) {
                ++j;
            }
            if (j == n) {
                std::copy(coords, coords + n, live.begin());
                extend(i);
                traced++;
                return;
            }
        }
        // the coordinates continuing the chain (c - 1 in |live|), by a merge: both are sorted
        scratch.clear();
        size_t k = 0;
        for (size_t j = 0; j < n; ++j) {
            const Coord c = coords[j];
            if (!c)
                continue;
            while (k < live.size() && live[k] < c - 1) {
                ++k;
            }
            if (k < live.size() && live[k] == c - 1)
                scratch.push_back(c);
        }
        if (scratch.empty()) {
            if (!live.empty())
                trace_breaks.push_back(i);  // supported, but the trace jumped
            // no coordinate continues the chain: start a new trace run here
            live.assign(coords, coords + n);
            start(i);
        } else {
            extend(i);
            live.assign(scratch.begin(), scratch.end());
        }
        traced++;
    }

    std::vector<KmerInterval> take_runs() {
        if (open.begin != kNone)
            runs.push_back(open);
        open = { kNone, kNone };
        return std::move(runs);
    }

  private:
    // the run ending at i continues through i; else a new run starts at i
    void extend(uint64_t i) {
        if (open.end == i) {
            open.end = i + 1;
        } else {
            start(i);
        }
    }
    void start(uint64_t i) {
        if (open.begin != kNone)
            runs.push_back(open);
        open = { i, i + 1 };
    }
};

// The size of the next batch of rows a /resolve decodes: the first kResolveFirstBatchRows
// rows, then about |target_bytes| of rows as wide as the widest of the batch before (rows
// widen and narrow along a query, conserved and variable regions, and the average would
// underestimate), at most twice that batch's rows (the next rows can be wider still) and at
// most |max_rows|
struct BatchSizer {
    size_t max_rows;
    uint64_t target_bytes;
    size_t next;

    BatchSizer(size_t max_rows, uint64_t target_bytes)
          : max_rows(std::max<size_t>(1, max_rows)), target_bytes(target_bytes),
            next(std::min(this->max_rows, kResolveFirstBatchRows)) {}

    template <class RowT>
    void done(const std::vector<RowT> &batch) {
        uint64_t widest = 1;
        for (const auto &row : batch) {
            widest = std::max(widest, annot::matrix::row_copy_bytes(row));
        }
        next = static_cast<size_t>(std::max<uint64_t>(1, std::min<uint64_t>({
                target_bytes / widest, 2 * batch.size(), max_rows })));
    }
};

std::vector<Row> anno_rows(const std::vector<node_index> &keys) {
    std::vector<Row> rows;
    rows.reserve(keys.size());
    for (node_index key : keys) {
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    return rows;
}

// The modelled bytes of an accumulator (SupportProfile::memory_peak): its fixed part, its runs
// and trace breaks (the chain's live coordinates are a handful, in place)
uint64_t support_bytes(const LabelSupport &s) {
    return sizeof(LabelSupport) + s.runs.size() * sizeof(KmerInterval)
         + s.trace_breaks.size() * sizeof(uint64_t);
}

template <class Supports>
uint64_t supports_bytes(const Supports &supports) {
    uint64_t bytes = 0;
    for (const LabelSupport &s : supports) {
        bytes += support_bytes(s);
    }
    return bytes;
}

// The heap of a chain's coordinates beyond those held in place (the discovery's charge)
uint64_t chain_heap_bytes(const LabelSupport &s) {
    const size_t capacity = s.live.capacity();
    return capacity > 4 ? capacity * sizeof(Coord) : 0;
}

uint64_t scratch_bytes(const std::vector<Coord> &scratch) {
    return scratch.capacity() * sizeof(Coord);
}

// What the budget has left (DecodeBudget::max_bytes() less what it holds)
uint64_t budget_left(const DecodeBudget &budget) {
    return budget.max_bytes() - budget.held();
}

// How a pass over the rows (for_each_row) ended
enum class PassEnd {
    COMPLETE,       // every key used
    CHECKPOINT,     // the checkpoint ended it between two batches (a deadline)
    MEMORY,         // the budget: a k-mer's row read alone, or its use, did not fit
};

// The budget's charge per key of the pass's setup (the next occurrence of each k-mer, kept for
// the pass) and of its transient hash table of the keys' last occurrences (seeded, 90 B a key:
// hopscotch_entry_bytes; T3 review R2, F1 and N1: 80 B, a table at most 90% full, was short of
// the 81.4 B a natural key set reached, and crafted node indices grew it to 604 B a key)
constexpr uint64_t kNextUseBytes = sizeof(size_t);
using KeyPositions = tsl::hopscotch_map<node_index, size_t, SeededHash>;
constexpr uint64_t kLastOccurrenceBytes
        = hopscotch_entry_bytes(sizeof(KeyPositions::value_type));
// a discovery's index of the records it tracks (header labels): the record, its accumulator's
// position (kDiscoverIndexEntryBytes an entry)
using RecordIndex = tsl::hopscotch_map<std::pair<Column, uint64_t>, uint64_t, SeededPairHash>;

template <class RowT>
uint64_t rows_bytes(const std::vector<RowT> &rows) {
    uint64_t bytes = 0;
    for (const RowT &row : rows) {
        bytes += annot::matrix::row_copy_bytes(row);
    }
    return bytes;
}

/**
 * Calls |use(j, row)| for every j of |keys| (the present k-mers' keys, in k-mer order) in
 * order, with the row of keys[j] decoded by |decode| in batches (BatchSizer: a batch is the
 * next rows to decode, not its k-mers). A key's row is decoded once per batch, and a row whose
 * key occurs again after its batch is kept for that occurrence while the rows kept take at
 * most |kept_max| bytes (row_copy_bytes), dropped at its last occurrence; beyond that bound a
 * repeated row is decoded again. |checkpoint| runs between batches; true from it ends the pass
 * there (a deadline). What the request holds is one batch of rows and the rows kept.
 *
 * Every byte the pass holds is charged to |budget| (SPEC-labeled-traversal-core.md §4.6): its
 * setup (the next occurrence of each key, and the transient table that finds them), each batch
 * by |decode(keys, &rows, &held)| — the rows of |keys| into |rows|, OK with the budget holding
 * |held[r]| for row r and IRowDiff::output_bytes<RowT>(keys.size()) for the vectors, or
 * REFUSED with nothing held and nothing returned — and each row kept with its entry (KeptRows:
 * a row's charge released at its last use, the entries' charge held as the most entries kept
 * at once while the table's buckets are, both released with the table when the rows are
 * dropped). |use(j, row)| returns false when the use of the row of keys[j] does not
 * fit beside what is held (nothing of it done; true: done, its own growth charged by it). A
 * charge that does not fit is met, in this order, by decoding the batch in halves, by giving
 * up the rows kept (decoded again when needed) and by reading the next row alone; the pass ends
 * (PassEnd::MEMORY) only at a k-mer whose row read alone, or whose use, does not fit beside what
 * the uses of the k-mers before it hold: that k-mer does not depend on the batches, the rows
 * kept or the time. With an unbounded budget nothing is refused and the pass is the one it was.
 *
 * Returns how many of |keys| were used: all of them, or those before the batch the pass ended
 * at (a checkpoint) or the k-mer it stopped at (the budget) — a prefix, since the keys are used
 * in order —, and sets |*end| to why. |rows_held_peak|, when given, is raised to the most row
 * bytes held at once (a batch and the rows kept, row_copy_bytes each:
 * SupportProfile::memory_peak)
 */
template <class RowT, class Decode, class Use, class Checkpoint>
size_t for_each_row(const std::vector<node_index> &keys, BatchSizer sizer, uint64_t kept_max,
                    DecodeBudget &budget, const Decode &decode, const Use &use,
                    const Checkpoint &checkpoint, PassEnd *end_cause,
                    uint64_t *rows_held_peak = nullptr) {
    constexpr size_t kNever = std::numeric_limits<size_t>::max();
    *end_cause = PassEnd::COMPLETE;
    // the next occurrence of each k-mer's key (kNever: none), held for the pass; the table of
    // the keys' last occurrences that builds it, held while it does
    const uint64_t next_use_bytes = annot::matrix::buffer_bytes(keys.size(), kNextUseBytes);
    const uint64_t setup_bytes = keys.size() * kLastOccurrenceBytes;
    if (!budget.charge(next_use_bytes + setup_bytes)) {
        *end_cause = PassEnd::MEMORY;
        return 0;
    }
    std::vector<size_t> next_use(keys.size(), kNever);
    {
        KeyPositions later;
        for (size_t j = keys.size(); j-- > 0; ) {
            auto [it, inserted] = later.try_emplace(keys[j], j);
            if (!inserted) {
                next_use[j] = it->second;
                it.value() = j;
            }
        }
    }
    budget.release(setup_bytes);
    // the rows kept for a later occurrence of their key (GPT review 6, #8: the table's buckets
    // outlive an entry's last use, so their charge is held until kept.drop() frees them)
    KeptRows<RowT, node_index> kept(budget);
    // the pass ends here, before keys[begin]: what it holds is released
    auto finish = [&](PassEnd cause) {
        *end_cause = cause;
        kept.drop();
        budget.release(next_use_bytes);
    };
    // the batch's keys (not charged, SPEC §4.6: at most 4,096 keys, 0.3 MiB at its growth, the
    // batch's vectors 0.16 MiB beside it; seeded, so that no key set grows it past that)
    KeyPositions in_batch;
    std::vector<node_index> batch_keys;
    std::vector<RowT> batch;
    std::vector<uint64_t> held;
    std::vector<bool> moved;
    for (size_t begin = 0; begin < keys.size(); ) {
        // the k-mers whose rows are kept or among the next |sizer.next| rows to decode
        in_batch.clear();
        batch_keys.clear();
        size_t end = begin;
        for ( ; end < keys.size(); ++end) {
            if (kept.contains(keys[end]) || in_batch.count(keys[end]))
                continue;
            if (batch_keys.size() == sizer.next)
                break;
            in_batch.emplace(keys[end], batch_keys.size());
            batch_keys.push_back(keys[end]);
        }
        batch.clear();
        held.clear();
        uint64_t containers = 0;
        if (!batch_keys.empty()) {
            if (decode(batch_keys, &batch, &held) != DecodeStatus::OK) {
                // the batch did not fit beside what is held: read it in halves; then without
                // the rows kept; a row that does not fit alone ends the pass at its k-mer (the
                // first of the batch: with no row kept, keys[begin] is batch_keys[0])
                if (batch_keys.size() > 1) {
                    sizer.next = batch_keys.size() / 2;
                } else if (kept.holds()) {
                    // (also a table whose rows were all used for the last time: its entries'
                    // charge is released with it, so that the stop does not depend on them)
                    kept.drop();
                } else {
                    finish(PassEnd::MEMORY);
                    return begin;
                }
                // a refused read took time too: the checkpoint is read before the next
                if (checkpoint()) {
                    finish(PassEnd::CHECKPOINT);
                    return begin;
                }
                continue;
            }
            containers = annot::matrix::IRowDiff::output_bytes<RowT>(batch_keys.size());
            sizer.done(batch);   // before any row is moved to |kept|
            // what the request holds now: this batch of rows and the rows kept
            if (rows_held_peak)
                *rows_held_peak = std::max(*rows_held_peak, kept.bytes() + rows_bytes(batch));
        }
        moved.assign(batch.size(), false);
        // the batch's rows not moved to |kept|, and its vectors
        auto release_batch = [&]() {
            for (size_t r = 0; r < batch.size(); ++r) {
                if (!moved[r])
                    budget.release(held[r]);
            }
            budget.release(containers);
        };
        size_t j = begin;
        for ( ; j < end; ++j) {
            auto kit = kept.find(keys[j]);
            if (kit != kept.end()) {
                if (!use(j, kit->second.row))
                    break;
                if (next_use[j] == kNever)
                    kept.last_use(kit);
                continue;
            }
            const size_t r = in_batch.find(keys[j])->second;
            RowT &row = batch[r];
            if (!use(j, row))
                break;
            // occurring again after this batch (none of its k-mers after j has the key)
            if (next_use[j] != kNever && next_use[j] >= end) {
                const uint64_t bytes = annot::matrix::row_copy_bytes(row);
                if (kept.bytes() + bytes <= kept_max
                        && kept.keep(keys[j], std::move(row), bytes, held[r])) {
                    moved[r] = true;
                }
            }
        }
        release_batch();
        if (j < end) {
            // the use of k-mer j did not fit: again from it without the rows kept, then with
            // its row alone; if neither leaves room, the pass ends before it
            begin = j;
            if (kept.holds()) {
                kept.drop();
            } else if (batch_keys.size() > 1) {
                sizer.next = 1;
            } else {
                finish(PassEnd::MEMORY);
                return begin;
            }
            if (checkpoint()) {
                finish(PassEnd::CHECKPOINT);
                return begin;
            }
            continue;
        }
        begin = end;
        // after the last batch the pass is complete whatever the checkpoint says
        if (checkpoint() && begin < keys.size()) {
            finish(PassEnd::CHECKPOINT);
            return begin;
        }
    }
    finish(PassEnd::COMPLETE);
    return keys.size();
}

} // namespace

// The refusals of resolve_support's options and query, before it reads anything
static void check_resolve_options(const LabelOracle &oracle,
                                  std::string_view query,
                                  const ResolveOptions &options) {
    if (options.labels.empty() == !options.discover) {
        throw std::invalid_argument("Specify either an explicit label list or discovery, not both");
    }
    // A discovery keeps at least one label, as the request parser requires: a cut keeping none
    // would have no last kept label to state its boundary by (LabelTruncation)
    if (options.discover && !options.discover_max_labels)
        throw std::invalid_argument("Discovery must keep at least one label");
    if (options.support == Support::TRACE) {
        if (!oracle.has_coordinates())
            throw std::invalid_argument("Trace support requires an annotation with k-mer coordinates");
        if (oracle.regime() != Regime::BASIC) {
            throw std::invalid_argument("Trace support is only defined for BASIC (forward-strand) "
                                        "graphs: coordinates carry no strand in canonical indexes");
        }
    }
    if (query.size() < oracle.get_k())
        throw std::invalid_argument("Query shorter than k");
}

void check_resolve_request(const LabelOracle &oracle,
                           std::string_view query,
                           const ResolveOptions &options) {
    check_resolve_options(oracle, query, options);
    // then, as resolve_support resolves them, the explicit labels or the header discovery's
    // record mapping
    if (!options.labels.empty()) {
        if (!options.drop_labels_not_in_graph) {
            for (const auto &name : options.labels) {
                oracle.resolve_label(name);
            }
        }
    } else if (options.discover_kind == LabelKind::HEADER && !oracle.coord_to_header()) {
        throw std::invalid_argument("Header discovery requires a CoordToHeader index");
    }
}

void check_explicit_seeds(const LabelOracle &oracle,
                          std::string_view query,
                          const ResolveOptions &options,
                          const SelectionPolicy &policy) {
    if (policy.policy != SelectionPolicy::EXPLICIT)
        return;
    const std::vector<ExplicitSeed> &seeds = policy.explicit_seeds;
    const uint64_t num_kmers = query.size() - oracle.get_k() + 1;
    auto out_of_range = [&](const ExplicitSeed &ex) {
        return ex.kmers.end > num_kmers || ex.kmers.begin >= ex.kmers.end;
    };
    // what a discovery profiles is known after the work, and a label dropped by a multi-graph
    // pair covers nothing rather than being refused: a seed's labels are then decided by the
    // work, and select_seeds checks them before any later seed
    const bool labels_after_work = options.discover || options.drop_labels_not_in_graph;
    // a seed's label that is no label of the request (explicit labels are all profiled): the
    // request's labels in one set, each seed's names looked up in it (a linear search per name
    // made a request of N labels and seeds naming them cost N^2 comparisons, review 5 #6)
    std::unordered_set<std::string_view> listed;
    if (!labels_after_work && !seeds.empty())
        listed.insert(options.labels.begin(), options.labels.end());
    auto unlisted = [&](const ExplicitSeed &ex) -> const std::string* {
        if (labels_after_work)
            return nullptr;
        for (const std::string &name : ex.labels) {
            if (!listed.count(name))
                return &name;
        }
        return nullptr;
    };
    size_t faulty = 0;
    while (faulty < seeds.size() && !out_of_range(seeds[faulty]) && !unlisted(seeds[faulty])) {
        // (the first seed whose labels the work decides ends what is certain here, after its
        // interval: a later seed's fault may not be select_seeds' first)
        if (labels_after_work && !seeds[faulty].labels.empty())
            return;
        ++faulty;
    }
    if (faulty == seeds.size())
        return;
    // A refusal for sure. Its words are select_seeds', whose first fault may be a seed not fully
    // in the graph before this one, or this one's before its labels: only such a request
    // maps the query's k-mers here (no row is read)
    const bool range = out_of_range(seeds[faulty]);
    if (faulty || !range) {
        const std::vector<node_index> keys = oracle.keys_of_sequence(query);
        const std::vector<KmerInterval> graph_runs
                = runs_of(keys.size(), [&](uint64_t i) { return keys[i] != npos; });
        for (size_t i = 0; i < faulty + !range; ++i) {
            const KmerInterval &seed = seeds[i].kmers;
            if (std::none_of(graph_runs.begin(), graph_runs.end(), [&](const KmerInterval &run) {
                    return run.begin <= seed.begin && run.end >= seed.end;
                })) {
                throw std::invalid_argument("Explicit seed interval is not fully present in the graph");
            }
        }
    }
    if (range)
        throw std::invalid_argument("Explicit seed interval out of range");
    throw std::invalid_argument(kSeedLabelNotProfiled
                                + ::utils::clip_for_message(*unlisted(seeds[faulty])) + "'");
}

SupportProfile resolve_support(LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();
    profile.support = options.support;

    check_resolve_options(oracle, query, options);

    // a phase boundary: the caller can abandon the request here (ResolveOptions::stop)
    auto checkpoint = [&]() {
        if (options.stop && options.stop() && options.abandon)
            options.abandon();
    };
    // A boundary of the work under a deadline (ResolveOptions::time_up): the client first, so
    // that a gone client is abandoned exactly where it was before; true ends the work here.
    // Without a deadline this is the checkpoint, at the same places
    auto poll = [&]() {
        checkpoint();
        return options.time_up && options.time_up();
    };
    // Where the deadline ended the work: the k-mers before |resolved| are resolved, the profile
    // being exactly the resolve of that prefix. The prefix is what makes a stopped answer honest
    // without a new statement per field: every run, count, truncation, candidate and seed is
    // that of a real query, the first |resolved| k-mers, never a sample of the whole query
    bool stopped = false;
    uint64_t resolved = 0;
    ResolveStop::Phase stop_phase = ResolveStop::ROWS;
    ResolveStop::Reason stop_reason = ResolveStop::TIME;
    bool cannot_progress = false;
    // the loops over the labels after the work: the answer's time, every kResolveCheckLabels
    auto finishing = [&](size_t i) {
        if (options.finish_check && i && i % kResolveCheckLabels == 0)
            options.finish_check();
    };
    profile.num_kmers = query.size() - profile.k + 1;
    checkpoint();
    std::vector<node_index> keys = oracle.keys_of_sequence(query);
    assert(keys.size() == profile.num_kmers);

    std::vector<bool> in_graph(profile.num_kmers);
    std::vector<node_index> present_keys;
    std::vector<uint64_t> present_pos;   // the k-mer index of each present key
    for (uint64_t i = 0; i < keys.size(); ++i) {
        in_graph[i] = keys[i] != npos;
        if (in_graph[i]) {
            present_keys.push_back(keys[i]);
            present_pos.push_back(i);
        }
    }
    profile.graph_runs = runs_of(in_graph);
    size_t num_present = present_keys.size();
    // The modelled memory (SupportProfile::memory_peak): what the query holds throughout — its
    // bases, every k-mer's key and presence, the present k-mers' keys and positions, the graph
    // runs — and what each phase adds to it at its largest
    const uint64_t query_held = query.size() + keys.size() * sizeof(node_index)
                              + (in_graph.size() + 7) / 8
                              + present_keys.size() * (sizeof(node_index) + sizeof(uint64_t))
                              + profile.graph_runs.size() * sizeof(KmerInterval);
    auto note_memory = [&](uint64_t bytes) {
        profile.memory_peak = std::max(profile.memory_peak, query_held + bytes);
    };
    note_memory(0);
    // the most the row passes held at once (a batch of rows and the rows kept)
    uint64_t rows_held = 0;
    // The work ended at k-mer x (|phase|): the profile becomes that of the query's first x
    // k-mers — num_kmers, and the graph runs cut there (the runs of the labels hold no k-mer
    // past x: only those before it were read). Applied once, before the labels' profiles are
    // taken, so that whatever is built from them is built from the prefix
    auto apply_stop = [&]() {
        if (!stopped || profile.stop)
            return;
        ResolveStop stop;
        stop.phase = stop_phase;
        stop.reason = stop_reason;
        stop.max_memory = stop_reason == ResolveStop::MEMORY ? options.discover_max_memory : 0;
        stop.cannot_progress = cannot_progress;
        stop.resolved_kmers = resolved;
        stop.query_kmers = profile.num_kmers;
        stop.query_graph_runs = profile.graph_runs;
        profile.stop = std::move(stop);
        profile.num_kmers = resolved;
        std::vector<KmerInterval> cut;
        for (const KmerInterval &run : profile.graph_runs) {
            if (run.begin >= resolved)
                break;
            cut.push_back({ run.begin, std::min(run.end, resolved) });
        }
        profile.graph_runs.swap(cut);
    };

    // ---- the rows: decoded in batches, each dropped once it is used (keeping them for a
    // second pass would hold about all of them, or decode them twice). |use(j, row)| gets the
    // row of present_keys[j], j ascending: tuple rows with |tuple_rows|, else whole rows
    // (for_each_row: a repeated key's row is kept for its later occurrences within
    // kept_bytes), and answers false when its accumulation does not fit |budget|. Returns how
    // many present k-mers were used: all, or under a deadline those before the batch it ended
    // at (none when it had passed before the first: |rows_expired|), or those before the k-mer
    // the memory bound stopped it at (|pass_end|)
    bool rows_expired = false;
    // A discovery's pass charges what it holds to one budget (ResolveOptions::
    // discover_max_memory, §4.6): unbounded without it, every read still the budget-aware one
    DecodeBudget budget(options.discover_max_memory ? options.discover_max_memory
                                                    : std::numeric_limits<uint64_t>::max());
    PassEnd pass_end = PassEnd::COMPLETE;
    // A batch's rows, charged to |budget|: by the budget-aware read where the annotation has
    // it (LabelOracle::decode_charged(): every buffer of the decode charged before it is
    // allocated, a read that does not fit refused whole), else by the default read, its rows
    // then charged as row_copy_bytes counts them (the batch can pass the bound while it is
    // decoded) and the read refused if they do not fit
    auto decode_rows = [&](const std::vector<node_index> &batch, auto *rows,
                           std::vector<uint64_t> *held) -> DecodeStatus {
        using RowT = typename std::decay_t<decltype(*rows)>::value_type;
        constexpr bool kTuples = std::is_same_v<RowT, MultiIntMatrix::RowTuples>;
        const std::vector<Row> anno = anno_rows(batch);
        std::vector<annot::matrix::RowCost> costs;
        DecodeStatus status;
        if constexpr (kTuples) {
            status = oracle.get_row_tuples(anno, budget, rows, &costs, held);
        } else {
            status = oracle.get_rows(anno, budget, rows, &costs, held);
        }
        if (status != DecodeStatus::UNSUPPORTED)
            return status;
        if constexpr (kTuples) {
            *rows = oracle.get_row_tuples(anno);
        } else {
            *rows = oracle.get_rows(anno);
        }
        held->resize(rows->size());
        uint64_t total = annot::matrix::IRowDiff::output_bytes<RowT>(rows->size());
        for (size_t r = 0; r < rows->size(); ++r) {
            (*held)[r] = annot::matrix::row_copy_bytes((*rows)[r]);
            total += (*held)[r];
        }
        if (!budget.charge(total)) {
            rows->clear();
            held->clear();
            return DecodeStatus::REFUSED;
        }
        return DecodeStatus::OK;
    };
    auto each_row = [&](bool tuple_rows, const auto &use_rows, const auto &use_tuples) -> size_t {
        if (rows_expired)
            return 0;
        const BatchSizer sizer(options.batch_rows, options.batch_bytes);
        if (tuple_rows) {
            return for_each_row<MultiIntMatrix::RowTuples>(present_keys, sizer, options.kept_bytes,
                budget, decode_rows, use_tuples, poll, &pass_end, &rows_held);
        } else {
            return for_each_row<BinaryMatrix::SetBitPositions>(present_keys, sizer,
                options.kept_bytes, budget, decode_rows, use_rows, poll, &pass_end, &rows_held);
        }
    };

    // ---- which labels to profile
    std::vector<LabelRef> refs;
    const bool with_coords = options.support == Support::TRACE;
    // a discovery's labels, with the support it accumulated while it read the rows
    std::vector<LabelSupport> discovered;
    if (!options.labels.empty()) {
        // each absent name listed once (a request may repeat a label)
        tsl::hopscotch_set<std::string_view> unknown;
        for (const auto &name : options.labels) {
            if (!options.drop_labels_not_in_graph) {
                refs.push_back(oracle.resolve_label(name));
            } else if (std::optional<LabelRef> ref = oracle.find_label(name)) {
                refs.push_back(std::move(*ref));
            } else if (unknown.insert(name).second) {
                profile.unknown_labels.push_back(name);
            }
        }
    } else {
        if (options.discover_kind == LabelKind::HEADER && !oracle.coord_to_header())
            throw std::invalid_argument("Header discovery requires a CoordToHeader index");
        // the deadline's first read: past it, no row is read
        rows_expired = poll();

        // Discovery reads every present k-mer's row once and accumulates, per label, both what
        // ranks it (its k-mers) and its profile (support runs; for trace the chain's state), so
        // that the top labels' profiles need no second read of the rows (a profile pass over
        // them made a column discover cost 1.6-2.1 times the profile of explicit labels on
        // refseq33m; keeping the rows for it held a query's worth of wide tuple rows). The
        // labels are those LabelQuery derives from the same rows: a column of the row (tuple
        // rows when trace needs its coordinates), and a sequence of a column holding one of
        // its coordinates (tuple rows), its local coordinates sorted.
        // The accumulators of the labels found, in the order they were first seen: a column's
        // by a slot per column id (4 bytes a column, where an accumulator per column of the
        // index, about 100 bytes, would be allocated and cleared by every request), a sequence's by
        // (column, seq_id) in one flat hash table (a std::map node per label took 0.8-1.7 us a
        // pair, 1.83 s for a header discover on 23S, refseq33m)
        const bool by_column = options.discover_kind == LabelKind::COLUMN;
        constexpr uint32_t kNoSlot = std::numeric_limits<uint32_t>::max();
        // The bounds (§4.6): what the labels' index and accumulators add is charged to
        // |budget| beside the rows — the slots once, per label tracked |label_bytes| (its
        // accumulator and id in deques, which never relocate, and a header label's entry in
        // the index), per run and trace break, and the heap of the chains' coordinates and of
        // the scratch as it grows. Before a k-mer's labels are accumulated, the most they can
        // add is checked against what is left (admit_row), so a charge after it always fits
        const bool bounded = options.discover_max_memory != 0;
        // a header discovery also counts the records it tracks per column, so that a k-mer's
        // new labels are bounded by those its columns have left (8 B a column)
        const uint64_t slot_bytes = annot::matrix::buffer_bytes(oracle.num_columns(),
                by_column ? sizeof(uint32_t) : sizeof(uint64_t));
        const bool index_fits = budget.charge(slot_bytes);
        std::vector<uint32_t> slot(by_column && index_fits ? oracle.num_columns() : 0, kNoSlot);
        std::vector<uint64_t> tracked_in(!by_column && index_fits ? oracle.num_columns() : 0, 0);
        // seeded (SeededPairHash): the records a query meets could otherwise be chosen to grow
        // the table past its charge (kDiscoverIndexEntryBytes; T3 review R2, F1)
        RecordIndex index;
        static_assert(kDiscoverIndexEntryBytes
                      == hopscotch_entry_bytes(sizeof(RecordIndex::value_type)));
        std::deque<LabelSupport> found;
        std::deque<std::pair<Column, uint64_t>> found_ids;
        const uint64_t label_bytes = sizeof(LabelSupport) + sizeof(std::pair<Column, uint64_t>)
                                   + sizeof(void*) + (by_column ? 0 : kDiscoverIndexEntryBytes);
        // The label cap (ResolveOptions::discover_max_tracked): once a k-mer's new labels would
        // pass it, no label first met there or later is admitted (|admitting| false); the
        // labels tracked keep being counted
        const size_t max_tracked = options.discover_max_tracked;
        bool admitting = true;
        LabelCap cap;
        size_t cut_index = 0;                // the cut k-mer's index among the present ones
        // a label's accumulator: created when it is met first while labels are admitted; null
        // for a label that is not tracked
        auto of = [&](Column c, uint64_t seq_id) -> LabelSupport* {
            if (by_column) {
                uint32_t &at = slot[c];
                if (at == kNoSlot) {
                    if (!admitting)
                        return nullptr;
                    at = static_cast<uint32_t>(found.size());
                    found.emplace_back();
                    found_ids.emplace_back(c, 0);
                }
                return &found[at];
            }
            if (!admitting) {
                auto it = index.find({ c, seq_id });
                return it != index.end() ? &found[it->second] : nullptr;
            }
            auto [it, inserted] = index.try_emplace({ c, seq_id }, found.size());
            if (inserted) {
                found.emplace_back();
                found_ids.emplace_back(c, seq_id);
                tracked_in[c]++;
            }
            return &found[it->second];
        };
        std::vector<Coord> scratch;
        const bool tuple_rows = with_coords || options.discover_kind == LabelKind::HEADER;
        std::vector<Coord> sorted;
        std::vector<Coord> local;
        // what the current k-mer's accumulation added to the runs, the trace breaks and the
        // chains' heap (charged after it, with its new labels and the scratch's growth)
        uint64_t grown = 0;
        auto add = [&](LabelSupport &s, uint64_t i) {
            const size_t runs = s.runs.size();
            s.add(i);
            grown += (s.runs.size() - runs) * kDiscoverRunBytes;
        };
        auto add_trace = [&](LabelSupport &s, uint64_t i, const Coord *coords, size_t n) {
            const size_t runs = s.runs.size();
            const size_t breaks = s.trace_breaks.size();
            const uint64_t heap = chain_heap_bytes(s);
            s.add_trace(i, coords, n, scratch);
            grown += (s.runs.size() - runs) * kDiscoverRunBytes
                   + (s.trace_breaks.size() - breaks) * kDiscoverBreakBytes
                   + (chain_heap_bytes(s) - heap);
        };
        // the scratch a trace chains and maps coordinates in (its growth charged after the k-mer,
        // within what admit_row checked)
        auto scratch_held = [&]() {
            return scratch_bytes(scratch) + scratch_bytes(local);
        };
        // The sort scratch, charged before it grows (its new buffer while the old one is still
        // held, the old one released once freed): false when the growth does not fit (nothing
        // changed). Only an unsorted column needs it (the row-diff forms' tuples are sorted)
        uint64_t sorted_charged = 0;
        auto reserve_sorted = [&](size_t n) {
            if (sorted.capacity() >= n)
                return true;
            if (!budget.charge(annot::matrix::buffer_bytes(n, sizeof(Coord))))
                return false;
            sorted.clear();
            sorted.reserve(n);
            budget.release(sorted_charged);
            sorted_charged = annot::matrix::buffer_bytes(n, sizeof(Coord));
            return true;
        };
        // a column's coordinates, sorted (a trace runs on sorted coordinates, as LabelQuery
        // sorts a hit's, and on sorted coordinates a sequence's are adjacent: the run mapper
        // maps a sequence's run for one rank/select); false when the sort scratch does not fit
        auto sorted_coords = [&](const auto &coords, const Coord **data) {
            if (std::is_sorted(coords.begin(), coords.end())) {
                *data = coords.data();
                return true;
            }
            if (!reserve_sorted(coords.size()))
                return false;
            sorted.assign(coords.begin(), coords.end());
            std::sort(sorted.begin(), sorted.end());
            *data = sorted.data();
            return true;
        };
        // The labels of a row not tracked yet and those tracked: lookups only (a column's
        // slot, a sequence's entry), each label once
        auto count_bits = [&](const BinaryMatrix::SetBitPositions &row) {
            uint64_t fresh = 0;
            for (Column c : row) {
                fresh += slot[c] == kNoSlot;
            }
            return std::make_pair(fresh, row.size() - fresh);
        };
        // (|sort_refused|: the sort scratch of the row's largest unsorted column, |unsorted|,
        // did not fit: nothing counted. It is grown to that column once, before any column is
        // sorted into it, so that no smaller growth is held beside it)
        bool sort_refused = false;
        auto count_tuples = [&](const MultiIntMatrix::RowTuples &row, size_t unsorted) {
            uint64_t fresh = 0;
            uint64_t known = 0;
            if (!by_column && !reserve_sorted(unsorted)) {
                sort_refused = true;
                return std::make_pair(fresh, known);
            }
            for (const auto &[c, coords] : row) {
                if (by_column) {
                    (slot[c] == kNoSlot ? fresh : known)++;
                    continue;
                }
                const Coord *data = nullptr;
                const bool sorted_fits = sorted_coords(coords, &data);  // (reserved above)
                assert(sorted_fits);
                (void)sorted_fits;
                LabelOracle::CoordRuns runs = oracle.coord_runs(c, coords.size());
                uint64_t previous = std::numeric_limits<uint64_t>::max();
                for (size_t t = 0; t < coords.size(); ++t) {
                    const uint64_t seq_id = runs.map(data[t]).first;
                    if (seq_id == previous)
                        continue;
                    previous = seq_id;
                    (index.count({ c, seq_id }) ? known : fresh)++;
                }
            }
            return std::make_pair(fresh, known);
        };
        // Whether the k-mer of present index j, whose row has at most |labels| labels, at most
        // |fresh_most| of them not tracked yet, can be accumulated: false when the most it can
        // add does not fit what |budget| has left (nothing done). |coord_bytes|: the most its
        // coordinates add (under trace the chains and their scratch; nothing under presence),
        // and |unsorted| its largest unsorted column, whose sort scratch it may grow. The labels
        // it would make tracked past the cap make it the cut k-mer (none of them admitted),
        // decided here and taken effect only once the k-mer fits. |count| gives the row's labels
        // not tracked and tracked (a pass of lookups, and for header labels the coordinates'
        // mapping), read only where the bounds alone do not decide
        auto admit_row = [&](size_t j, uint64_t labels, uint64_t fresh_most, uint64_t coord_bytes,
                             size_t unsorted, const auto &count) {
            // (a flag and a value rather than an optional, whose copy GCC 13 takes for
            // maybe-uninitialized)
            bool is_counted = false;
            std::pair<uint64_t, uint64_t> counted { 0, 0 };
            auto counts = [&]() {
                if (!is_counted) {
                    counted = count();
                    is_counted = true;
                }
                return counted;
            };
            sort_refused = false;
            bool cut_here = false;
            uint64_t fresh_max = 0;
            if (admitting) {
                fresh_max = fresh_most;
                if (max_tracked && found.size() + fresh_most > max_tracked) {
                    fresh_max = counts().first;
                    if (sort_refused)
                        return false;
                    if (found.size() + fresh_max > max_tracked) {
                        cut_here = true;
                        fresh_max = 0;
                    }
                }
            }
            if (bounded) {
                // the new labels, a run and a break per label touched, the coordinates, and the
                // sort scratch the accumulation grows (charged as it grows)
                auto most = [&](uint64_t fresh, uint64_t touched) {
                    const uint64_t sort = sorted.capacity() < unsorted
                            ? annot::matrix::buffer_bytes(unsorted, sizeof(Coord)) : 0;
                    return fresh * label_bytes
                         + touched * (kDiscoverRunBytes + kDiscoverBreakBytes)
                         + coord_bytes + sort;
                };
                if (most(fresh_max, std::min<uint64_t>(labels, found.size() + fresh_max))
                        > budget_left(budget)) {
                    const auto [fresh, known] = counts();
                    if (sort_refused)
                        return false;
                    const uint64_t admitted = (cut_here || !admitting) ? 0 : fresh;
                    if (most(admitted, known + admitted) > budget_left(budget))
                        return false;
                }
            }
            if (cut_here) {
                admitting = false;
                cut_index = j;
                cap.cut_kmer = present_pos[j];
                cap.new_at_cut = counts().first;
                cap.labels_at_cut = counts().first + counts().second;
                cap.tracked = found.size();
            }
            return true;
        };
        // after a k-mer's accumulation: its new labels and growth, within what admit_row saw
        // left (so the charge fits)
        auto charge_row = [&](size_t labels_before, uint64_t scratch_before) {
            const uint64_t bytes = (found.size() - labels_before) * label_bytes + grown
                                 + (scratch_held() - scratch_before);
            grown = 0;
            const bool charged = budget.charge(bytes);
            assert(charged && "within the bound admit_row checked");
            (void)charged;
        };
        const size_t used = !index_fits ? 0 : each_row(tuple_rows,
            [&](size_t j, const BinaryMatrix::SetBitPositions &row) {
                // a column not tracked yet is one of those the index has left
                const uint64_t fresh_most = std::min<uint64_t>(
                        row.size(), oracle.num_columns() - found.size());
                if (!admit_row(j, row.size(), fresh_most, 0, 0, [&]() { return count_bits(row); }))
                    return false;
                const size_t labels_before = found.size();
                const uint64_t scratch_before = scratch_held();
                for (Column c : row) {
                    if (LabelSupport *s = of(c, 0))
                        add(*s, present_pos[j]);
                }
                charge_row(labels_before, scratch_before);
                return true;
            },
            [&](size_t j, const MultiIntMatrix::RowTuples &row) {
                // a column label per column, a sequence label per coordinate at most, and of the
                // labels not tracked yet at most those the index has left: per column (header
                // labels: its records not tracked yet), or of all columns
                const annot::CoordToHeader *cth = oracle.coord_to_header();
                uint64_t coords_in_row = 0;
                uint64_t fresh_most = 0;
                size_t widest = 0;
                size_t unsorted = 0;
                for (const auto &[c, coords] : row) {
                    coords_in_row += coords.size();
                    widest = std::max(widest, coords.size());
                    if (!std::is_sorted(coords.begin(), coords.end()))
                        unsorted = std::max(unsorted, coords.size());
                    if (!by_column) {
                        // (a column the record mapping does not have: its coordinates bound it,
                        // and the mapping refuses them as the base did)
                        fresh_most += c < cth->num_columns()
                                ? std::min<uint64_t>(coords.size(),
                                                     cth->num_sequences(c) - tracked_in[c])
                                : coords.size();
                    }
                }
                if (by_column)
                    fresh_most = std::min<uint64_t>(row.size(), oracle.num_columns() - found.size());
                // what the coordinates add: under trace each one's place in a chain, and the
                // scratch they are chained in (and for header labels mapped in) grown to twice
                // the widest column; under presence nothing (the sort scratch apart)
                const uint64_t scratch_grows = std::max<uint64_t>(
                        2 * widest, scratch.capacity()) - scratch.capacity()
                        + (by_column ? 0 : std::max<uint64_t>(2 * widest, local.capacity())
                                               - local.capacity());
                const uint64_t coord_bytes = with_coords
                        ? coords_in_row * kDiscoverChainBytes + scratch_grows * sizeof(Coord)
                        : 0;
                const uint64_t scratch_before = scratch_held();
                if (!admit_row(j, by_column ? row.size() : coords_in_row, fresh_most, coord_bytes,
                               unsorted, [&]() { return count_tuples(row, unsorted); }))
                    return false;
                // the sort scratch at the row's largest unsorted column, once (within what
                // admit_row checked), so that every column below sorts within it
                const bool sort_fits = reserve_sorted(unsorted);
                assert(sort_fits && "within the bound admit_row checked");
                (void)sort_fits;
                const size_t labels_before = found.size();
                const uint64_t i = present_pos[j];
                for (const auto &[c, coords] : row) {
                    const Coord *data = nullptr;
                    const bool sorted_fits = sorted_coords(coords, &data);  // (reserved above)
                    assert(sorted_fits);
                    (void)sorted_fits;
                    if (options.discover_kind == LabelKind::COLUMN) {
                        // tuple rows only for a trace (presence reads whole rows)
                        if (LabelSupport *s = of(c, 0))
                            add_trace(*s, i, data, coords.size());
                        continue;
                    }
                    uint64_t previous = std::numeric_limits<uint64_t>::max();
                    if (!with_coords) {
                        // each sequence of the column counted once per k-mer
                        oracle.map_coords(c, data, coords.size(),
                                          [&, column = c](Coord, uint64_t seq_id, Coord) {
                            if (seq_id != previous) {
                                if (LabelSupport *s = of(column, seq_id))
                                    add(*s, i);
                            }
                            previous = seq_id;
                        });
                        continue;
                    }
                    // each sequence's local coordinates, in order (sorted)
                    local.clear();
                    auto trace = [&](Column column, uint64_t seq_id) {
                        if (LabelSupport *s = of(column, seq_id))
                            add_trace(*s, i, local.data(), local.size());
                    };
                    oracle.map_coords(c, data, coords.size(),
                                      [&, column = c](Coord, uint64_t seq_id, Coord at) {
                        if (seq_id != previous && !local.empty()) {
                            trace(column, previous);
                            local.clear();
                        }
                        previous = seq_id;
                        local.push_back(at);
                    });
                    if (!local.empty())
                        trace(c, previous);
                }
                charge_row(labels_before, scratch_before);
                return true;
            });
        if (!index_fits)
            pass_end = PassEnd::MEMORY;
        // the pass at its largest: the next occurrence of every present k-mer, the rows held,
        // and the accumulators (taken at their end, where they are largest) with the labels'
        // slots or ids
        note_memory(present_keys.size() * sizeof(size_t) + rows_held + supports_bytes(found)
                    + found_ids.size() * sizeof(found_ids[0])
                    + (by_column ? slot.size() * sizeof(uint32_t)
                                 : index.size() * (sizeof(std::pair<Column, uint64_t>)
                                                   + 2 * sizeof(uint64_t))));
        if (used < num_present) {
            // The deadline or the memory bound ended the pass: the rows of the first |used|
            // present k-mers were read, so the k-mers before the next present one are resolved
            // (those between are absent from the graph, which the mapping told). The labels met
            // are those of that prefix, each with its support in it: what a discovery on the
            // prefix accumulates (under the cap: the labels first met before the cut k-mer, a
            // k-mer of the prefix)
            stopped = true;
            stop_phase = ResolveStop::ROWS;
            stop_reason = pass_end == PassEnd::MEMORY ? ResolveStop::MEMORY : ResolveStop::TIME;
            // the memory bound at the query's first k-mer in the graph: the sequence from there
            // holds the same k-mers in the graph, so the same setup and the same first row, and
            // a resolve of it stops there again under the same bound
            cannot_progress = stop_reason == ResolveStop::MEMORY && used == 0;
            resolved = present_pos[used];
            num_present = used;
        }
        checkpoint();
        // A label found, as the ranking orders it: its k-mers supported, where its accumulator
        // and id are (found, found_ids), and its name once it is looked up (null until then).
        // The k-mers supported are the number the answer states as kmers_supported (SPEC-
        // multigraph-envelope.md §7): under support "trace" the traced k-mers (those with a
        // coordinate of the label, add_trace), else the k-mers it is present at. Both are
        // final here: the pass over the rows is over (whole, or the prefix a deadline left,
        // num_present cut to it), and nothing after it adds to an accumulator's counts
        struct Ranked {
            uint64_t kmers;
            uint64_t at;
            const std::string *name;
        };
        std::vector<Ranked> ranked;
        ranked.reserve(found.size());
        for (uint64_t at = 0; at < found.size(); ++at) {
            finishing(at);
            ranked.push_back({ with_coords ? found[at].traced : found[at].kmers, at, nullptr });
        }
        std::vector<uint32_t>().swap(slot);
        RecordIndex().swap(index);
        // More k-mers supported first, ties by the name the answer states (bytes: the column's name, or
        // the sequence's header), then by (column, seq_id) for a header that two columns hold.
        // A tie by column id let the index's numbering of its columns decide which of the
        // tied labels a cut keeps; by name the same labels rank the same in every index that
        // holds them, so that a client can merge the discoveries of several graphs by
        // (support, name), and each cut states the tie group it went through
        // (LabelTruncation::tied_at_boundary, SPEC-multigraph-envelope.md §7).
        std::sort(ranked.begin(), ranked.end(),
                  [](const Ranked &a, const Ranked &b) { return a.kmers > b.kmers; });
        const size_t kept = std::min(ranked.size(), options.discover_max_labels);
        // The labels named: the kept ones and the rest of the last kept label's tie group.
        // A name is a reference into the label encoder or the headers (no decoding), but
        // ordering every label found by name would read a name per comparison, scattered
        // over the headers, where the k-mer counts are adjacent integers — and a header
        // discovery can find far more sequences than it keeps. So only the labels the order
        // needs are named: the kept ones, which the answer names anyway (the order of their
        // equal-support runs is the order it states), and the dropped part of the boundary's
        // tie group, to choose which of the group are kept. The labels dropped below the
        // boundary's support are never named.
        size_t named = kept;
        while (named < ranked.size() && ranked[named].kmers == ranked[kept - 1].kmers) {
            named++;
        }
        for (size_t i = 0; i < named; ++i) {
            finishing(i);
            const auto &[column, seq_id] = found_ids[ranked[i].at];
            ranked[i].name = by_column ? &oracle.column_name(column)
                                       : &oracle.header_name(column, seq_id);
        }
        auto by_name = [&](const Ranked &a, const Ranked &b) {
            const int order = a.name->compare(*b.name);
            return order ? order < 0 : found_ids[a.at] < found_ids[b.at];
        };
        // each run of equal support in name order; the boundary's run only as far as it is
        // kept (its dropped part needs no order)
        size_t boundary = 0;
        for (size_t begin = 0; begin < kept; ) {
            size_t end = begin + 1;
            while (end < named && ranked[end].kmers == ranked[begin].kmers) {
                end++;
            }
            if (end <= kept) {
                std::sort(ranked.begin() + begin, ranked.begin() + end, by_name);
            } else {
                std::partial_sort(ranked.begin() + begin, ranked.begin() + kept,
                                  ranked.begin() + end, by_name);
            }
            boundary = begin;
            begin = end;
        }
        if (ranked.size() > kept) {
            LabelTruncation trunc;
            trunc.total = ranked.size();
            trunc.kept = kept;
            trunc.min_kept_kmers = ranked[kept - 1].kmers;
            trunc.max_dropped_kmers = ranked[kept].kmers;
            // full length: supported at every k-mer of the query present in the graph (of
            // the prefix read, under a stop) -- under trace traced at every one of them
            for (size_t i = kept; i < ranked.size(); ++i) {
                if (ranked[i].kmers == num_present)
                    trunc.dropped_full_length++;
            }
            trunc.tied_at_boundary = named - boundary;
            trunc.tied_kept = kept - boundary;
            profile.labels_truncated = trunc;
        }
        if (!admitting) {
            // The cap bound: the labels tracked are those first met before the cut k-mer, each
            // counted over the whole query (of the prefix, under a stop); an untracked label is
            // present only at the in-graph k-mers from the cut k-mer on, so it supports at most
            // that many. The kept list is the uncut discovery's top max_labels when that many
            // were kept and the last of them supports more (an untracked label cannot then
            // reach it, nor tie with it)
            cap.max_tracked = max_tracked;
            assert(cap.tracked == found.size());
            cap.untracked_max_kmers = num_present - cut_index;
            cap.max_labels = options.discover_max_labels;
            cap.kept = kept;
            cap.min_kept_kmers = kept ? ranked[kept - 1].kmers : 0;
            cap.ranking_exact = kept == options.discover_max_labels
                                && cap.min_kept_kmers > cap.untracked_max_kmers;
            profile.labels_capped = cap;
        }
        ranked.resize(kept);
        discovered.reserve(ranked.size());
        for (const Ranked &label : ranked) {
            finishing(refs.size());
            LabelRef ref;
            ref.kind = options.discover_kind;
            ref.column = found_ids[label.at].first;
            ref.seq_id = found_ids[label.at].second;
            ref.name = *label.name;
            refs.push_back(std::move(ref));
            discovered.push_back(std::move(found[label.at]));
        }
    }

    profile.labels.reserve(refs.size());
    for (const auto &ref : refs) {
        LabelProfile lp;
        lp.label = ref;
        profile.labels.push_back(lp);
    }
    // a discovery that met no label (or whose deadline passed before its first row)
    apply_stop();
    if (refs.empty())
        return profile;

    if (options.labels.empty()) {
        // a discovery's profiles are what it accumulated
        for (LabelId l = 0; l < refs.size(); ++l) {
            finishing(l);
            LabelProfile &lp = profile.labels[l];
            LabelSupport &support = discovered[l];
            lp.runs = support.take_runs();
            lp.trace_breaks = std::move(support.trace_breaks);
            lp.kmers_supported = with_coords ? support.traced : support.kmers;
        }
    } else {
        // ---- support per k-mer of explicit labels
        LabelQuery query_labels(oracle, refs, with_coords);
        // Every label's support from one pass over the k-mers: each hit goes to its label's
        // accumulator, k-mer by k-mer in ascending order, as a discovery accumulates them. Under
        // trace a per-label scan (for each label, every k-mer's hit list searched for it) would
        // cost O(labels x k-mers x hits), an absent label reading every list whole: 2,000
        // explicit labels took 5.9 s that way where the discovery returning the same profiles
        // takes 0.2 s, with nothing polling the client meanwhile. Each label receives the calls
        // such a scan makes, in the same order, so the profiles are the same. A presence bitmap
        // of labels x k-mers would take 125 MB for 1,000 labels on 1 M k-mers; the accumulator's
        // runs are runs_of's of that bitmap
        std::vector<LabelSupport> support(refs.size());
        std::vector<Coord> scratch;
        auto scatter = [&](uint64_t i, const LabelQuery::NodeHits &node_hits) {
            for (const auto &h : node_hits) {
                assert(h.label < refs.size());
                LabelSupport &s = support[h.label];
                if (with_coords) {
                    // trace-consistent runs: a chain of coordinates increasing by one per
                    // k-mer. Column labels have coordinates in the column frame, header labels
                    // in the sequence frame; either way consecutive k-mers must have
                    // consecutive coords (LabelSupport::add_trace). A label's first hit of the
                    // k-mer is its one (the scan stopped at it; a query holds a label once)
                    if (s.last == i)
                        continue;
                    s.add_trace(i, h.coords.data(), h.coords.size(), scratch);
                } else {
                    // present at the k-mer once, however many hits name it (as its bit was)
                    if (s.open.end == i + 1)
                        continue;
                    s.add(i);
                }
            }
        };
        // The row paths under a deadline: the hits of the k-mers [scattered, end), whose keys
        // the priming has all cached (no row is read: a fetch decodes only missing keys), in
        // pieces of kResolveCheckKmers k-mers with the client checked between them
        uint64_t scattered = 0;
        std::vector<node_index> piece;
        auto scatter_primed = [&](uint64_t end) {
            for (uint64_t begin = scattered; begin < end; begin += kResolveCheckKmers) {
                if (begin > scattered)
                    checkpoint();
                const uint64_t e = std::min<uint64_t>(end, begin + kResolveCheckKmers);
                piece.assign(keys.begin() + begin, keys.begin() + e);
                const auto hits = query_labels.fetch(piece);
                for (uint64_t i = begin; i < e; ++i) {
                    scatter(i, hits[i - begin]);
                }
            }
            scattered = std::max(scattered, end);
        };
        // The deadline's first read: the labels are resolved (an unknown one refused), no row
        // is read yet. |resolved| is the end of the k-mers whose hits the query can give
        // without reading a row it has not read: all of them until a stop
        const bool expired = poll();
        resolved = profile.num_kmers;
        const std::string path = query_labels.access_path();
        if (expired) {
            // what the rows (or a direct read's cells) would have been read for: k-mers
            // before the first present one need none
            stopped = num_present > 0;
            stop_phase = path != "direct" ? ResolveStop::ROWS : ResolveStop::SUPPORT;
            resolved = num_present ? present_pos[0] : profile.num_kmers;
        }
        // The distinct keys' whole or tuple rows are decoded in batches and primed into the
        // query (the hits a fetch builds from them, kept per key), so that the fetch below
        // decodes nothing and the request holds one batch of rows at a time beside the hits —
        // each row once, as the fetch decoded the distinct keys; direct cell reads hold no rows
        if (path != "direct" && !expired) {
            std::vector<node_index> distinct;
            // under a deadline, the k-mer where each distinct key first occurs: the keys are
            // primed in that order, so a stop after the first P of them has primed every k-mer
            // before distinct_at[P] (a key occurring before it occurred first before it)
            std::vector<uint64_t> distinct_at;
            {
                tsl::hopscotch_set<node_index, SeededHash> seen;
                for (size_t j = 0; j < present_keys.size(); ++j) {
                    if (seen.insert(present_keys[j]).second) {
                        distinct.push_back(present_keys[j]);
                        if (options.time_up)
                            distinct_at.push_back(present_pos[j]);
                    }
                }
            }
            BatchSizer sizer(options.batch_rows, options.batch_bytes);
            std::vector<node_index> batch;
            // the distinct keys (and under a deadline where each occurs first) beside the query
            const uint64_t distinct_held = distinct.size() * sizeof(node_index)
                                         + distinct_at.size() * sizeof(uint64_t);
            for (size_t begin = 0; begin < distinct.size(); begin += batch.size()) {
                batch.assign(distinct.begin() + begin,
                             distinct.begin() + std::min(distinct.size(), begin + sizer.next));
                // a batch of rows held while its hits are primed into the query
                if (path == "tuples") {
                    const auto rows = oracle.get_row_tuples(anno_rows(batch));
                    query_labels.prime(batch, rows);
                    sizer.done(rows);
                    note_memory(distinct_held + rows_bytes(rows) + query_labels.cache_bytes()
                                + supports_bytes(support));
                } else {
                    const auto rows = oracle.get_rows(anno_rows(batch));
                    query_labels.prime(batch, rows);
                    sizer.done(rows);
                    note_memory(distinct_held + rows_bytes(rows) + query_labels.cache_bytes()
                                + supports_bytes(support));
                }
                const size_t next = begin + batch.size();
                if (options.time_up) {
                    // every k-mer before distinct_at[next] has its key primed now (the whole
                    // query after the last batch): their hits are scattered here, so that a
                    // stop between batches keeps every row it read (a separate hits pass
                    // would read the deadline again and cut the answer short of the rows
                    // already primed)
                    scatter_primed(next < distinct.size() ? distinct_at[next] : resolved);
                }
                // after the last batch the rows are all read, whatever the deadline says
                if (poll() && next < distinct.size()) {
                    stopped = true;
                    stop_phase = ResolveStop::ROWS;
                    resolved = distinct_at[next];
                    break;
                }
            }
        }

        if (!options.time_up) {
            // no deadline: one fetch of every k-mer's hits
            auto hits = query_labels.fetch(keys);
            // every k-mer's hits held at once beside the query's cache of them
            uint64_t hits_held = hits.size() * sizeof(hits[0]);
            for (const auto &h : hits) {
                hits_held += LabelQuery::held_bytes(h);
            }
            note_memory(hits_held + query_labels.cache_bytes());
            checkpoint();
            for (uint64_t i = 0; i < hits.size(); ++i) {
                // a client gone is seen within a few thousand k-mers of the pass, not after it
                if (i && i % kResolveCheckKmers == 0)
                    checkpoint();
                scatter(i, hits[i]);
            }
        } else if (path != "direct") {
            // the row paths scattered their hits with the priming. Left, when the deadline
            // passed before the first row: the k-mers before the first present one, absent
            // from the graph, with no hits
            scatter_primed(resolved);
        } else {
            // Under a deadline the direct path's hits are fetched kResolveCheckKmers k-mers at a
            // time, the deadline read before each piece but the first (read just before it): one
            // fetch of the whole query is a piece no clock read can end — it reads every k-mer's
            // cell of every label. A key's hits do not depend on the piece it is fetched in, so
            // the profile is the one fetch's
            for (uint64_t begin = 0; begin < resolved; begin += kResolveCheckKmers) {
                if (begin && poll()) {
                    // resolved: the k-mers before the first present one from |begin| on (any
                    // between are absent); none before |resolved|, nothing was left to read
                    const auto next = std::lower_bound(present_pos.begin(), present_pos.end(),
                                                       begin);
                    if (next != present_pos.end() && *next < resolved) {
                        stopped = true;
                        stop_phase = ResolveStop::SUPPORT;
                        resolved = *next;
                    }
                    break;
                }
                const uint64_t end = std::min<uint64_t>(resolved, begin + kResolveCheckKmers);
                piece.assign(keys.begin() + begin, keys.begin() + end);
                const auto hits = query_labels.fetch(piece);
                for (uint64_t i = begin; i < end; ++i) {
                    scatter(i, hits[i - begin]);
                }
            }
        }
        apply_stop();
        // the accumulators at their end, beside the query's cache of hits
        note_memory(supports_bytes(support) + query_labels.cache_bytes());
        for (LabelId l = 0; l < refs.size(); ++l) {
            finishing(l);
            LabelProfile &lp = profile.labels[l];
            lp.runs = support[l].take_runs();
            // presence keeps none (add records no trace)
            lp.trace_breaks = std::move(support[l].trace_breaks);
            lp.kmers_supported = with_coords ? support[l].traced : support[l].kmers;
        }
    }

    // ---- seed candidates: identical maximal runs grouped
    std::map<KmerInterval, std::vector<LabelId>> groups;
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        finishing(l);
        for (const auto &run : profile.labels[l].runs) {
            if (run.size() >= options.min_block_kmers)
                groups[run].push_back(l);
        }
    }
    for (auto &[iv, labels] : groups) {
        std::sort(labels.begin(), labels.end());
        profile.candidates.push_back({ iv, std::move(labels) });
    }
    std::sort(profile.candidates.begin(), profile.candidates.end(),
              [](const SeedCandidate &a, const SeedCandidate &b) {
                  if (a.kmers.size() != b.kmers.size())
                      return a.kmers.size() > b.kmers.size();
                  if (a.labels.size() != b.labels.size())
                      return a.labels.size() > b.labels.size();
                  return a.kmers.begin < b.kmers.begin;
              });
    // the profile, complete, beside the query
    note_memory(profile_bytes(profile));
    return profile;
}

uint64_t profile_bytes(const SupportProfile &p) {
    uint64_t bytes = sizeof(SupportProfile) + p.graph_runs.size() * sizeof(KmerInterval);
    for (const LabelProfile &lp : p.labels) {
        bytes += sizeof(LabelProfile) + lp.label.name.size()
               + lp.runs.size() * sizeof(KmerInterval) + lp.trace_breaks.size() * sizeof(uint64_t);
    }
    for (const SeedCandidate &c : p.candidates) {
        bytes += sizeof(SeedCandidate) + c.labels.size() * sizeof(LabelId);
    }
    for (const std::string &name : p.unknown_labels) {
        bytes += sizeof(std::string) + name.size();
    }
    if (p.stop)
        bytes += p.stop->query_graph_runs.size() * sizeof(KmerInterval);
    return bytes;
}

uint64_t selection_bytes(const SeedSelection &selection) {
    auto names = [](const std::vector<std::string> &list) {
        uint64_t bytes = list.size() * sizeof(std::string);
        for (const std::string &name : list) {
            bytes += name.size();
        }
        return bytes;
    };
    uint64_t bytes = sizeof(SeedSelection);
    for (const FrozenSeed &seed : selection.seeds) {
        bytes += sizeof(FrozenSeed) + seed.seed_id.size() + seed.sequence.size()
               + names(seed.labels) + names(seed.population.dropped)
               + seed.population.dropped_digest.size()
               + seed.overlaps_with.size() * sizeof(size_t) + names(seed.labels_not_covering);
    }
    return bytes;
}


/********************************* selection *********************************/

namespace {

// A Fenwick tree of counts over positions 0..n-1: add, the sum over positions < i, and
// the position of the k-th unit, each in O(log n).
struct Fenwick {
    std::vector<uint32_t> tree;
    size_t n;
    size_t top = 1;   // the largest power of two <= n (1 when n == 0)

    explicit Fenwick(size_t size) : tree(size + 1, 0), n(size) {
        while (top * 2 <= n) {
            top *= 2;
        }
    }
    void add(size_t i, uint32_t v) {
        for (++i; i <= n; i += i & (~i + 1)) {
            tree[i] += v;
        }
    }
    uint32_t prefix(size_t i) const {   // sum over positions < i
        uint32_t s = 0;
        for (; i > 0; i -= i & (~i + 1)) {
            s += tree[i];
        }
        return s;
    }
    size_t kth(uint32_t k) const {   // the 0-based position holding the k-th unit (1-based k)
        size_t pos = 0;
        for (size_t step = top; step; step >>= 1) {
            if (pos + step <= n && tree[pos + step] < k) {
                pos += step;
                k -= tree[pos];
            }
        }
        return pos;
    }
};

// does label l have a run containing [begin, end)?
bool covers(const LabelProfile &lp, const KmerInterval &iv) {
    for (const auto &run : lp.runs) {
        if (run.begin <= iv.begin && run.end >= iv.end)
            return true;
    }
    return false;
}

std::vector<LabelId> covering_labels(const SupportProfile &profile, const KmerInterval &iv) {
    std::vector<LabelId> out;
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        if (covers(profile.labels[l], iv))
            out.push_back(l);
    }
    return out;
}

struct Picked {
    KmerInterval kmers;
    std::vector<LabelId> labels;
    std::vector<std::string> labels_not_covering;
    size_t rank = 0;  // the order in which the policy picked it
};

} // namespace

SeedSelection select_seeds(const SupportProfile &profile,
                           std::string_view query,
                           const SelectionPolicy &policy_in,
                           bool canonical_orientation) {
    SeedSelection selection;
    selection.policy = policy_in;
    SelectionPolicy &policy = selection.policy;
    if (policy.min_block_bp < profile.k)
        policy.min_block_bp = profile.k;
    const uint64_t min_kmers = policy.min_block_bp - profile.k + 1;
    selection.num_candidates = profile.candidates.size();

    std::vector<Picked> picked;

    switch (policy.policy) {
        case SelectionPolicy::LONGEST_FIRST: {
            for (const auto &cand : profile.candidates) {
                if (cand.kmers.size() < min_kmers)
                    continue;
                selection.num_eligible++;
                if (picked.size() < policy.max_seeds)
                    picked.push_back({ cand.kmers, cand.labels, {} });
            }
            break;
        }
        case SelectionPolicy::MAX_SUPPORT: {
            // Maximize the number of labels whose support run *contains* the chosen
            // interval [a, b): a run [s, e) qualifies iff s <= a and e >= b.
            //
            // For a fixed a, the count |{e >= b}| is non-increasing in b, so the best
            // count is reached at the shortest admissible b = a + min_kmers, and among
            // the intervals with that count the longest one ends at the smallest
            // qualifying end, i.e. b = min{e : e >= a + min_kmers}. So one binary
            // search per candidate a suffices. Both conditions only change at
            // endpoints, so a ranges over run starts and over (e - min_kmers).
            // This is O(R log R) per round for R runs, instead of scanning pairs.
            struct Run { uint64_t begin, end; };
            std::vector<Run> runs;
            for (const auto &lp : profile.labels) {
                for (const auto &run : lp.runs) {
                    if (run.size() >= min_kmers)
                        runs.push_back({ run.begin, run.end });
                }
            }
            const std::vector<Run> all_runs = [&runs]() {
                std::sort(runs.begin(), runs.end(),
                          [](const Run &a, const Run &b) { return a.begin < b.begin; });
                return runs;
            }();

            // intervals already claimed by a picked seed, as sorted disjoint ranges
            std::vector<KmerInterval> taken;

            while (picked.size() < policy.max_seeds) {
                // Rebuild the runs with the claimed intervals subtracted, so later rounds
                // see clipped runs and their endpoints. Reusing the original endpoints and
                // merely rejecting overlapping intervals loses the best remaining interval
                // (e.g. with A=[0,10), B=[0,5) and [0,5) taken, [5,10) is only reachable
                // through the clipped run's own endpoints).
                runs.clear();
                for (const Run &run : all_runs) {
                    uint64_t begin = run.begin;
                    for (const auto &t : taken) {
                        if (t.end <= begin || t.begin >= run.end)
                            continue;
                        if (t.begin > begin && t.begin - begin >= min_kmers)
                            runs.push_back({ begin, t.begin });
                        begin = std::max(begin, t.end);
                    }
                    if (begin < run.end && run.end - begin >= min_kmers)
                        runs.push_back({ begin, run.end });
                }
                std::sort(runs.begin(), runs.end(),
                          [](const Run &a, const Run &b) { return a.begin < b.begin; });
                if (runs.empty())
                    break;

                std::vector<uint64_t> candidate_a;
                for (const auto &run : runs) {
                    candidate_a.push_back(run.begin);
                    if (run.end >= min_kmers)
                        candidate_a.push_back(run.end - min_kmers);
                }
                std::sort(candidate_a.begin(), candidate_a.end());
                candidate_a.erase(std::unique(candidate_a.begin(), candidate_a.end()),
                                  candidate_a.end());

                std::optional<Picked> best;
                size_t best_count = 0;
                size_t next_run = 0;
                // the ends of the runs with begin <= a, as counts over the compressed run
                // ends: an insertion, "how many ends >= x" and "the smallest end >= x" each
                // cost O(log r) (keeping them in a sorted vector made the sweep O(r^2))
                std::vector<uint64_t> coords;
                coords.reserve(runs.size());
                for (const Run &run : runs) {
                    coords.push_back(run.end);
                }
                std::sort(coords.begin(), coords.end());
                coords.erase(std::unique(coords.begin(), coords.end()), coords.end());
                auto index_of = [&](uint64_t x) {
                    return static_cast<size_t>(std::lower_bound(coords.begin(), coords.end(), x)
                                               - coords.begin());
                };
                Fenwick ends(coords.size());
                uint32_t inserted = 0;
                for (uint64_t a : candidate_a) {
                    while (next_run < runs.size() && runs[next_run].begin <= a) {
                        ends.add(index_of(runs[next_run].end), 1);
                        inserted++;
                        next_run++;
                    }
                    if (a + min_kmers > profile.num_kmers)
                        break;
                    const size_t lo = index_of(a + min_kmers);
                    const uint32_t below = lo < coords.size() ? ends.prefix(lo) : inserted;
                    if (below == inserted)
                        continue;
                    KmerInterval iv { a, coords[ends.kth(below + 1)] };
                    size_t count = inserted - below;
                    bool better = !best
                        || count > best_count
                        || (count == best_count
                            && (iv.size() > best->kmers.size()
                                || (iv.size() == best->kmers.size() && iv.begin < best->kmers.begin)));
                    if (better) {
                        // only the count is needed to compare candidates; the label
                        // list is materialised once, for the winner, below (doing it on
                        // every improvement made a round O(r^2) in the number of runs)
                        best = Picked{ iv, {}, {} };
                        best_count = count;
                    }
                }
                if (!best)
                    break;
                best->labels = covering_labels(profile, best->kmers);
                assert(best->labels.size() == best_count);
                taken.push_back(best->kmers);
                std::sort(taken.begin(), taken.end());
                picked.push_back(*best);
            }
            selection.num_eligible = picked.size();
            break;
        }
        case SelectionPolicy::EXPLICIT: {
            // The seeds' names looked up in one map of the profiled labels (a name's FIRST
            // label, as a scan in profile order finds it) and one set of the labels the pair
            // does not hold: a linear search per name made seeds naming N labels of a profile of
            // N cost N^2 comparisons (review 5 #6)
            std::unordered_map<std::string_view, size_t> profiled;
            profiled.reserve(profile.labels.size());
            for (size_t l = 0; l < profile.labels.size(); ++l) {
                profiled.emplace(profile.labels[l].label.name, l);
            }
            const std::unordered_set<std::string_view> absent(profile.unknown_labels.begin(),
                                                              profile.unknown_labels.end());
            for (const auto &ex : policy.explicit_seeds) {
                if (ex.kmers.end > profile.num_kmers || ex.kmers.begin >= ex.kmers.end)
                    throw std::invalid_argument("Explicit seed interval out of range");
                bool in_one_run = false;
                for (const auto &run : profile.graph_runs) {
                    if (run.begin <= ex.kmers.begin && run.end >= ex.kmers.end)
                        in_one_run = true;
                }
                if (!in_one_run)
                    throw std::invalid_argument("Explicit seed interval is not fully present in the graph");
                Picked p;
                p.kmers = ex.kmers;
                for (const auto &name : ex.labels) {
                    const auto it = profiled.find(name);
                    if (it == profiled.end()) {
                        // a label this pair does not hold (unknown_labels, a multi-graph
                        // server's) covers nothing of the seed here
                        if (absent.count(name)) {
                            p.labels_not_covering.push_back(name);
                            continue;
                        }
                        throw std::invalid_argument(kSeedLabelNotProfiled
                                                    + ::utils::clip_for_message(name) + "'");
                    }
                    if (covers(profile.labels[it->second], ex.kmers)) {
                        p.labels.push_back(it->second);
                    } else {
                        p.labels_not_covering.push_back(name);
                    }
                }
                selection.num_eligible++;
                if (!p.labels.empty())
                    picked.push_back(std::move(p));
            }
            break;
        }
    }

    // Merge overlapping picks onto their common interval, carried by the labels of both.
    // Only LONGEST_FIRST picks can overlap: MAX_SUPPORT picks outside the earlier ones.
    // A merged pick keeps the rank of its best member, and the picks stay in rank order.
    if (policy.merge_overlapping && policy.policy == SelectionPolicy::LONGEST_FIRST) {
        for (size_t i = 0; i < picked.size(); ++i) {
            picked[i].rank = i;
        }
        std::sort(picked.begin(), picked.end(),
                  [](const Picked &a, const Picked &b) { return a.kmers < b.kmers; });
        std::vector<Picked> merged;
        for (auto &p : picked) {
            if (!merged.empty() && merged.back().kmers.end > p.kmers.begin) {
                Picked &m = merged.back();
                KmerInterval common { std::max(m.kmers.begin, p.kmers.begin),
                                      std::min(m.kmers.end, p.kmers.end) };
                if (common.size() >= min_kmers) {
                    m.kmers = common;
                    m.labels = covering_labels(profile, common);
                    m.rank = std::min(m.rank, p.rank);
                    continue;
                }
            }
            merged.push_back(std::move(p));
        }
        std::sort(merged.begin(), merged.end(),
                  [](const Picked &a, const Picked &b) { return a.rank < b.rank; });
        picked.swap(merged);
    }

    // freeze
    for (auto &p : picked) {
        FrozenSeed seed;
        seed.kmers = p.kmers;
        seed.sequence = std::string(query.substr(p.kmers.begin, p.kmers.size() + profile.k - 1));
        seed.labels_not_covering = p.labels_not_covering;

        std::vector<std::pair<std::string, LabelId>> ordered;
        for (LabelId l : p.labels) {
            ordered.emplace_back(profile.labels[l].label.name, l);
        }
        switch (policy.label_order) {
            case SelectionPolicy::HASH:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    std::string sa = a.first + "\t" + std::to_string(policy.sample_seed);
                    std::string sb = b.first + "\t" + std::to_string(policy.sample_seed);
                    auto ha = fnv1a64(sa), hb = fnv1a64(sb);
                    return std::tie(ha, a.first) < std::tie(hb, b.first);
                });
                break;
            case SelectionPolicy::COLUMN_ID:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    const auto &ra = profile.labels[a.second].label, &rb = profile.labels[b.second].label;
                    return std::tie(ra.column, ra.seq_id) < std::tie(rb.column, rb.seq_id);
                });
                break;
            case SelectionPolicy::KMERS_SUPPORTED:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    auto ka = profile.labels[a.second].kmers_supported;
                    auto kb = profile.labels[b.second].kmers_supported;
                    return std::tie(kb, a.first) < std::tie(ka, b.first);  // more first
                });
                break;
        }
        seed.population.supporting_total = ordered.size();
        seed.population.order = policy.label_order;
        seed.population.sample_seed = policy.sample_seed;
        uint64_t digest = 0xcbf29ce484222325ULL;
        for (size_t i = 0; i < ordered.size(); ++i) {
            if (i < policy.max_labels_per_seed) {
                seed.labels.push_back(ordered[i].first);
            } else {
                seed.population.dropped_count++;
                digest = fnv1a64(ordered[i].first + "\n", digest);
                if (ordered.size() - policy.max_labels_per_seed <= 1000)
                    seed.population.dropped.push_back(ordered[i].first);
            }
        }
        seed.population.included = seed.labels.size();
        if (seed.population.dropped_count)
            seed.population.dropped_digest = hex64(digest);
        std::sort(seed.labels.begin(), seed.labels.end());
        seed.seed_id = make_seed_id(policy.release_id, seed.sequence, canonical_orientation, seed.labels);
        selection.seeds.push_back(std::move(seed));
    }

    for (size_t i = 0; i < selection.seeds.size(); ++i) {
        for (size_t j = 0; j < selection.seeds.size(); ++j) {
            if (i == j)
                continue;
            const auto &a = selection.seeds[i].kmers, &b = selection.seeds[j].kmers;
            if (a.begin < b.end && b.begin < a.end)
                selection.seeds[i].overlaps_with.push_back(j);
        }
    }
    return selection;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
