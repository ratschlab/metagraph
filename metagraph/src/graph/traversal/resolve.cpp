#include "resolve.hpp"

#include <algorithm>
#include <iterator>
#include <cctype>
#include <map>
#include <numeric>
#include <tuple>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "graph/annotated_dbg.hpp"
#include "annotation/binary_matrix/row_diff/row_diff_cache.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;
using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;

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

/**
 * Calls |use(j, row)| for every j of |keys| (the present k-mers' keys, in k-mer order) in
 * order, with the row of keys[j] decoded by |decode| (a vector of keys to their rows) in
 * batches (BatchSizer: a batch is the next rows to decode, not its k-mers). A key's row is
 * decoded once per batch, and a row whose key occurs again after its batch is kept for that
 * occurrence while the rows kept take at most |kept_max| bytes (row_copy_bytes), dropped at its
 * last occurrence; beyond that bound a repeated row is decoded again. |checkpoint| runs between
 * batches; true from it ends the pass there (a deadline). What the request holds is one batch
 * of rows and the rows kept. Returns how many of |keys| were used: all of them, or those before
 * the batch the pass ended at — a prefix, since the keys are used in order
 */
template <class RowT, class Decode, class Use, class Checkpoint>
size_t for_each_row(const std::vector<node_index> &keys, BatchSizer sizer, uint64_t kept_max,
                    const Decode &decode, const Use &use, const Checkpoint &checkpoint) {
    constexpr size_t kNever = std::numeric_limits<size_t>::max();
    // the next occurrence of each k-mer's key (kNever: none)
    std::vector<size_t> next_use(keys.size(), kNever);
    {
        tsl::hopscotch_map<node_index, size_t> later;
        for (size_t j = keys.size(); j-- > 0; ) {
            auto [it, inserted] = later.try_emplace(keys[j], j);
            if (!inserted) {
                next_use[j] = it->second;
                it.value() = j;
            }
        }
    }
    tsl::hopscotch_map<node_index, std::pair<RowT, uint64_t>> kept;
    uint64_t kept_bytes = 0;
    tsl::hopscotch_map<node_index, size_t> in_batch;
    std::vector<node_index> batch_keys;
    for (size_t begin = 0; begin < keys.size(); ) {
        // the k-mers whose rows are kept or among the next |sizer.next| rows to decode
        in_batch.clear();
        batch_keys.clear();
        size_t end = begin;
        for ( ; end < keys.size(); ++end) {
            if (kept.count(keys[end]) || in_batch.count(keys[end]))
                continue;
            if (batch_keys.size() == sizer.next)
                break;
            in_batch.emplace(keys[end], batch_keys.size());
            batch_keys.push_back(keys[end]);
        }
        std::vector<RowT> batch;
        if (!batch_keys.empty()) {
            batch = decode(batch_keys);
            sizer.done(batch);   // before any row is moved to |kept|
        }
        for (size_t j = begin; j < end; ++j) {
            auto kit = kept.find(keys[j]);
            if (kit != kept.end()) {
                use(j, kit->second.first);
                if (next_use[j] == kNever) {
                    kept_bytes -= kit->second.second;
                    kept.erase(kit);
                }
                continue;
            }
            RowT &row = batch[in_batch.find(keys[j])->second];
            use(j, row);
            // occurring again after this batch (none of its k-mers after j has the key)
            if (next_use[j] != kNever && next_use[j] >= end) {
                const uint64_t bytes = annot::matrix::row_copy_bytes(row);
                if (kept_bytes + bytes <= kept_max) {
                    kept_bytes += bytes;
                    kept.emplace(keys[j], std::make_pair(std::move(row), bytes));
                }
            }
        }
        begin = end;
        // after the last batch the pass is complete whatever the checkpoint says
        if (checkpoint() && begin < keys.size())
            return begin;
    }
    return keys.size();
}

} // namespace

SupportProfile resolve_support(LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();
    profile.support = options.support;

    if (options.labels.empty() == !options.discover) {
        throw std::invalid_argument("Specify either an explicit label list or discovery, not both");
    }
    if (options.support == Support::TRACE) {
        if (!oracle.has_coordinates())
            throw std::invalid_argument("Trace support requires an annotation with k-mer coordinates");
        if (oracle.regime() != Regime::BASIC) {
            throw std::invalid_argument("Trace support is only defined for BASIC (forward-strand) "
                                        "graphs: coordinates carry no strand in canonical indexes");
        }
    }
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");

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
    // The work ended at k-mer x (|phase|): the profile becomes that of the query's first x
    // k-mers — num_kmers, and the graph runs cut there (the runs of the labels hold no k-mer
    // past x: only those before it were read). Applied once, before the labels' profiles are
    // taken, so that whatever is built from them is built from the prefix
    auto apply_stop = [&]() {
        if (!stopped || profile.stop)
            return;
        profile.stop = ResolveStop { stop_phase, resolved, profile.num_kmers,
                                     profile.graph_runs };
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
    // kept_bytes). Returns how many present k-mers were used: all, or under a deadline those
    // before the batch it ended at (none when it had passed before the first: |rows_expired|)
    bool rows_expired = false;
    auto each_row = [&](bool tuple_rows, const auto &use_rows, const auto &use_tuples) -> size_t {
        if (rows_expired)
            return 0;
        const BatchSizer sizer(options.batch_rows, options.batch_bytes);
        if (tuple_rows) {
            return for_each_row<MultiIntMatrix::RowTuples>(present_keys, sizer, options.kept_bytes,
                [&](const std::vector<node_index> &batch) {
                    return oracle.get_row_tuples(anno_rows(batch));
                }, use_tuples, poll);
        } else {
            return for_each_row<BinaryMatrix::SetBitPositions>(present_keys, sizer, options.kept_bytes,
                [&](const std::vector<node_index> &batch) {
                    return oracle.get_rows(anno_rows(batch));
                }, use_rows, poll);
        }
    };

    // ---- which labels to profile
    std::vector<LabelRef> refs;
    const bool with_coords = options.support == Support::TRACE;
    // a discovery's labels, with the support it accumulated while it read the rows
    std::vector<LabelSupport> discovered;
    if (!options.labels.empty()) {
        for (const auto &name : options.labels) {
            refs.push_back(oracle.resolve_label(name));
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
        std::vector<uint32_t> slot(by_column ? oracle.num_columns() : 0, kNoSlot);
        tsl::hopscotch_map<std::pair<Column, uint64_t>, uint64_t, LabelKeyHash> index;
        std::vector<LabelSupport> found;
        std::vector<std::pair<Column, uint64_t>> found_ids;
        auto of = [&](Column c, uint64_t seq_id) -> LabelSupport& {
            if (by_column) {
                uint32_t &at = slot[c];
                if (at == kNoSlot) {
                    at = static_cast<uint32_t>(found.size());
                    found.emplace_back();
                    found_ids.emplace_back(c, 0);
                }
                return found[at];
            }
            auto [it, inserted] = index.try_emplace({ c, seq_id }, found.size());
            if (inserted) {
                found.emplace_back();
                found_ids.emplace_back(c, seq_id);
            }
            return found[it->second];
        };
        std::vector<Coord> scratch;
        const bool tuple_rows = with_coords || options.discover_kind == LabelKind::HEADER;
        std::vector<Coord> sorted;
        std::vector<Coord> local;
        const size_t used = each_row(tuple_rows,
            [&](size_t j, const BinaryMatrix::SetBitPositions &row) {
                for (Column c : row) {
                    of(c, 0).add(present_pos[j]);
                }
            },
            [&](size_t j, const MultiIntMatrix::RowTuples &row) {
                const uint64_t i = present_pos[j];
                for (const auto &[c, coords] : row) {
                    // a trace runs on sorted coordinates (LabelQuery sorts a hit's), and on
                    // sorted coordinates a sequence's are adjacent: the run mapper
                    // (LabelOracle::map_coords) maps a sequence's run for one rank/select
                    const Coord *data = coords.data();
                    if (!std::is_sorted(coords.begin(), coords.end())) {
                        sorted.assign(coords.begin(), coords.end());
                        std::sort(sorted.begin(), sorted.end());
                        data = sorted.data();
                    }
                    if (options.discover_kind == LabelKind::COLUMN) {
                        // tuple rows only for a trace (presence reads whole rows)
                        of(c, 0).add_trace(i, data, coords.size(), scratch);
                        continue;
                    }
                    uint64_t previous = std::numeric_limits<uint64_t>::max();
                    if (!with_coords) {
                        // each sequence of the column counted once per k-mer
                        oracle.map_coords(c, data, coords.size(),
                                          [&, column = c](Coord, uint64_t seq_id, Coord) {
                            if (seq_id != previous)
                                of(column, seq_id).add(i);
                            previous = seq_id;
                        });
                        continue;
                    }
                    // each sequence's local coordinates, in order (sorted)
                    local.clear();
                    oracle.map_coords(c, data, coords.size(),
                                      [&, column = c](Coord, uint64_t seq_id, Coord at) {
                        if (seq_id != previous && !local.empty()) {
                            of(column, previous).add_trace(i, local.data(), local.size(),
                                                           scratch);
                            local.clear();
                        }
                        previous = seq_id;
                        local.push_back(at);
                    });
                    if (!local.empty())
                        of(c, previous).add_trace(i, local.data(), local.size(), scratch);
                }
            });
        if (used < num_present) {
            // The deadline ended the pass: the rows of the first |used| present k-mers were
            // read, so the k-mers before the next present one are resolved (those between are
            // absent from the graph, which the mapping told). The labels met are those of that
            // prefix, each with its support in it: what a discovery on the prefix accumulates
            stopped = true;
            stop_phase = ResolveStop::ROWS;
            resolved = present_pos[used];
            num_present = used;
        }
        checkpoint();
        std::vector<std::tuple<std::pair<Column, uint64_t>, uint64_t, uint64_t>> ranked;
        ranked.reserve(found.size());
        for (uint64_t at = 0; at < found.size(); ++at) {
            finishing(at);
            ranked.emplace_back(found_ids[at], found[at].kmers, at);
        }
        std::vector<uint32_t>().swap(slot);
        tsl::hopscotch_map<std::pair<Column, uint64_t>, uint64_t, LabelKeyHash>().swap(index);
        // more k-mers first; ties by column id, then seq_id (the order the std::map gave)
        std::sort(ranked.begin(), ranked.end(), [](const auto &a, const auto &b) {
            return std::get<1>(a) != std::get<1>(b) ? std::get<1>(a) > std::get<1>(b)
                                                    : std::get<0>(a) < std::get<0>(b);
        });
        if (ranked.size() > options.discover_max_labels) {
            LabelTruncation trunc;
            trunc.total = ranked.size();
            trunc.kept = options.discover_max_labels;
            trunc.min_kept_kmers = std::get<1>(ranked[trunc.kept - 1]);
            trunc.max_dropped_kmers = std::get<1>(ranked[trunc.kept]);
            for (size_t i = trunc.kept; i < ranked.size(); ++i) {
                if (std::get<1>(ranked[i]) == num_present)
                    trunc.dropped_full_length++;
            }
            ranked.resize(trunc.kept);
            profile.labels_truncated = trunc;
        }
        discovered.reserve(ranked.size());
        for (const auto &[id, count, at] : ranked) {
            finishing(refs.size());
            LabelRef ref;
            ref.kind = options.discover_kind;
            ref.column = id.first;
            ref.seq_id = id.second;
            ref.name = ref.kind == LabelKind::COLUMN ? oracle.column_name(ref.column)
                                                     : oracle.header_name(ref.column, ref.seq_id);
            refs.push_back(ref);
            discovered.push_back(std::move(found[at]));
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
                tsl::hopscotch_set<node_index> seen;
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
            for (size_t begin = 0; begin < distinct.size(); begin += batch.size()) {
                batch.assign(distinct.begin() + begin,
                             distinct.begin() + std::min(distinct.size(), begin + sizer.next));
                if (path == "tuples") {
                    const auto rows = oracle.get_row_tuples(anno_rows(batch));
                    query_labels.prime(batch, rows);
                    sizer.done(rows);
                } else {
                    const auto rows = oracle.get_rows(anno_rows(batch));
                    query_labels.prime(batch, rows);
                    sizer.done(rows);
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
    return profile;
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
                    auto it = std::find_if(profile.labels.begin(), profile.labels.end(),
                                           [&](const LabelProfile &lp) { return lp.label.name == name; });
                    if (it == profile.labels.end())
                        throw std::invalid_argument("Explicit seed label was not profiled: '" + name + "'");
                    if (covers(*it, ex.kmers)) {
                        p.labels.push_back(it - profile.labels.begin());
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

    // merge overlapping picks onto their common interval
    if (policy.merge_overlapping && policy.policy != SelectionPolicy::EXPLICIT) {
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
                    continue;
                }
            }
            merged.push_back(std::move(p));
        }
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
