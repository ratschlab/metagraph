#include "label_oracle.hpp"

#include <algorithm>
#include <cmath>

#include <tsl/hopscotch_set.h>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;
using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;
using annot::matrix::DecodeBudget;
using annot::matrix::DecodeStatus;
using annot::matrix::RowCost;
using annot::matrix::buffer_bytes;
using annot::matrix::small_vector_bytes;

// What a budget-aware cache entry is taken to hold besides its value's own buffers: its share
// of the hash map's buckets (a slot of key and value at load 0.4 after a growth) with the
// transient of a rehash (both bucket arrays), and its entry in the costs map likewise. The
// caches stay within their byte bound in these terms (an estimate of the map's slack, the
// values' buffers as decode_budget.hpp models them).
constexpr uint64_t kCacheEntryBytes = 256;
constexpr uint64_t kCostEntryBytes = 128;
// the work of a key's dependency rows: as a fetched row's (8 per row, 1 per entry)
static uint64_t dependency_units(const RowCost &cost) {
    return 8 * cost.dependency_rows + cost.dependency_entries;
}


LabelOracle::LabelOracle(const AnnotatedDBG &anno_graph,
                         const annot::CoordToHeader *coord_to_header_override)
      : anno_graph_(anno_graph) {
    const DeBruijnGraph &loaded = anno_graph_.get_graph();

    // Regime detection: a CanonicalDBG wrapper means a PRIMARY index. Otherwise the
    // mode tells native CANONICAL graphs from BASIC ones.
    if (const auto *wrapper = dynamic_cast<const CanonicalDBG*>(&loaded)) {
        regime_ = Regime::PRIMARY;
        // Per-request clone: same node ids (offset = primary max_index), own caches.
        // The shared wrapper is never mutated.
        local_canonical_ = std::make_shared<CanonicalDBG>(wrapper->get_graph_ptr());
        dbg_succ_ = dynamic_cast<const DBGSuccinct*>(&wrapper->get_graph());
        if (dbg_succ_) {
            node_first_cache_ = std::make_shared<NodeFirstCache>(*dbg_succ_);
            local_canonical_->add_extension(node_first_cache_);
        }
        canonical_ = local_canonical_.get();
        graph_ = canonical_;
    } else {
        if (loaded.get_mode() == DeBruijnGraph::PRIMARY) {
            throw std::invalid_argument("PRIMARY graphs must be wrapped into CanonicalDBG "
                                        "before traversal");
        }
        regime_ = loaded.get_mode() == DeBruijnGraph::CANONICAL ? Regime::CANONICAL
                                                               : Regime::BASIC;
        graph_ = &loaded;
        dbg_succ_ = dynamic_cast<const DBGSuccinct*>(&loaded);
        if (dbg_succ_)
            node_first_cache_ = std::make_shared<NodeFirstCache>(*dbg_succ_);
        sshash_ = dynamic_cast<const DBGSSHash*>(&loaded);
    }

    // get_matrix() may lock and flush (ColumnCompressed): call it once per request.
    matrix_ = &anno_graph_.get_annotator().get_matrix();
    get_entry_ = dynamic_cast<const annot::matrix::GetEntrySupport*>(matrix_);
    tuples_ = dynamic_cast<const MultiIntMatrix*>(matrix_);
    rd_ = dynamic_cast<const annot::matrix::IRowDiff*>(matrix_);
    coord_to_header_ = coord_to_header_override ? coord_to_header_override
                                                : anno_graph_.get_coord_to_header();
    if (coord_to_header_ && !tuples_) {
        logger->warn("CoordToHeader present but the annotation has no coordinates; "
                     "header labels are unavailable");
        coord_to_header_ = nullptr;
    }
}

LabelOracle::~LabelOracle() {}

size_t LabelOracle::get_k() const { return graph_->get_k(); }

uint64_t LabelOracle::num_columns() const { return matrix_->num_columns(); }
uint64_t LabelOracle::num_rows() const { return matrix_->num_rows(); }

std::vector<node_index> LabelOracle::keys_of_sequence(std::string_view sequence) const {
    std::vector<node_index> keys;
    if (sequence.size() < get_k())
        return keys;
    keys.reserve(sequence.size() - get_k() + 1);
    // map_to_nodes returns annotation keys in every regime: the node itself (BASIC),
    // the primary node (CanonicalDBG), the canonical node (native CANONICAL graphs)
    graph_->map_to_nodes(sequence, [&](node_index key) { keys.push_back(key); });
    counters_.keys_mapped += keys.size();
    return keys;
}

node_index LabelOracle::key_of(node_index node, std::string_view kmer) const {
    counters_.keys_mapped++;
    if (node == npos)
        return npos;
    switch (regime_) {
        case Regime::BASIC:
            return node;
        case Regime::PRIMARY:
            return canonical_->get_base_node(node);
        case Regime::CANONICAL:
            if (sshash_)
                return std::min(node, static_cast<const DBGSSHash&>(*graph_).reverse_complement(node));
            assert(kmer.size() == get_k());
            {
                node_index key = npos;
                graph_->map_to_nodes(kmer, [&](node_index k) { key = k; });
                return key;
            }
    }
    return npos;
}

std::vector<node_index>
LabelOracle::keys_of_path(const std::vector<node_index> &nodes, std::string_view window) const {
    assert(window.size() == nodes.size() + get_k() - 1);
    std::vector<node_index> keys(nodes.size());
    switch (regime_) {
        case Regime::BASIC:
            keys = nodes;
            break;
        case Regime::PRIMARY:
            for (size_t i = 0; i < nodes.size(); ++i) {
                keys[i] = nodes[i] == npos ? npos : canonical_->get_base_node(nodes[i]);
            }
            break;
        case Regime::CANONICAL:
            if (sshash_) {
                const auto &sshash = static_cast<const DBGSSHash&>(*graph_);
                for (size_t i = 0; i < nodes.size(); ++i) {
                    keys[i] = nodes[i] == npos ? npos
                                               : std::min(nodes[i], sshash.reverse_complement(nodes[i]));
                }
            } else {
                // one incremental mapping over the whole window instead of O(k) per node
                size_t i = 0;
                graph_->map_to_nodes(window, [&](node_index key) { keys[i++] = key; });
                assert(i == nodes.size());
            }
            break;
    }
    counters_.keys_mapped += nodes.size();
    return keys;
}

std::optional<Column> LabelOracle::find_column(const std::string &name) const {
    const auto &encoder = anno_graph_.get_annotator().get_label_encoder();
    if (!encoder.label_exists(name))
        return std::nullopt;
    return static_cast<Column>(encoder.encode(name));
}

const std::string& LabelOracle::column_name(Column column) const {
    return anno_graph_.get_annotator().get_label_encoder().decode(column);
}

std::optional<LabelRef> LabelOracle::find_header(const std::string &name) const {
    if (!coord_to_header_)
        return std::nullopt;
    // the reverse index belongs to the CoordToHeader (built once per loaded index, and
    // gone with it), never to a cache that could outlive it
    auto found = coord_to_header_->find_header(name);
    if (!found)
        return std::nullopt;
    LabelRef ref;
    ref.kind = LabelKind::HEADER;
    ref.column = found->first;
    ref.seq_id = found->second;
    ref.name = name;
    return ref;
}

const std::string& LabelOracle::header_name(Column column, uint64_t seq_id) const {
    assert(coord_to_header_);
    return coord_to_header_->get_headers(column).at(seq_id);
}

std::pair<uint64_t, Coord> LabelOracle::map_coord(Column column, Coord coord) const {
    assert(coord_to_header_);
    counters_.coords_mapped++;
    return coord_to_header_->map_single_coord(column, coord);
}

uint64_t LabelOracle::num_kmers_in_sequence(Column column, uint64_t seq_id) const {
    assert(coord_to_header_);
    return coord_to_header_->num_kmers_in_sequence(column, seq_id);
}

LabelRef LabelOracle::resolve_label(const std::string &name) const {
    if (auto column = find_column(name)) {
        LabelRef ref;
        ref.kind = LabelKind::COLUMN;
        ref.column = *column;
        ref.name = name;
        return ref;
    }
    if (auto header = find_header(name))
        return *header;
    throw std::invalid_argument("Label not found in the annotation: '" + name + "'"
                                + (coord_to_header_ ? "" : " (no sequence header index loaded)"));
}

size_t DecodePacer::next(size_t remaining, double ms_left, size_t previous,
                         double previous_ms) const {
    if (!(target_ms > 0) || !remaining)
        return remaining;
    // the per-row time to plan with: this read's own previous chunk once there is one, else
    // the slowest the request has seen (0: none measured yet, or too fast to measure)
    const bool known = previous > 0 || ms_per_row > 0;
    const double rate = previous ? previous_ms / static_cast<double>(previous) : ms_per_row;
    const double factor = previous ? rest_factor : far_factor;
    // The deadline cannot fall into the rest of the read: one piece, as before pass 5 (the
    // products are NaN for a rate of 0 with an infinite factor, which then splits). Nothing
    // measured yet: one piece only without a deadline, else a first chunk measures the rows
    if (known ? static_cast<double>(remaining) * rate * factor < ms_left
              : std::isinf(ms_left) && !std::isinf(factor))
        return remaining;
    const double budget = std::min(target_ms, ms_left);
    double n;
    if (!(budget > 0)) {
        n = 1;      // the deadline has passed or is now: the caller's stop decides
    } else if (rate > 0) {
        n = budget / rate;
    } else {
        n = known ? std::numeric_limits<double>::infinity() : static_cast<double>(first_rows);
    }
    // every split read starts small and grows by at most 4 times per chunk, so that rows slower
    // than any measured before are met by a small chunk
    n = std::min(n, previous ? 4.0 * static_cast<double>(previous)
                             : static_cast<double>(first_rows));
    if (!(n >= 1))
        return 1;
    return n >= static_cast<double>(remaining) ? remaining : static_cast<size_t>(n);
}

void DecodePacer::record(size_t rows, double ms) {
    max_read_ms = std::max(max_read_ms, ms);
    if (!rows || !(ms >= 0))
        return;
    // The slowest per-row time seen, kept for the request: it only decides whether a read is
    // one piece or starts with a small chunk, so a pessimistic value costs a first chunk at
    // most, while an optimistic one lets a whole read overrun the deadline. Row costs vary by
    // two orders of magnitude on a row-diff annotation (a row whose path is shared with the
    // other rows of its call against one that decodes its own): a rate that decayed to the
    // cheap rows let a call of 182 unshared rows run 287 ms past a 200 ms budget (review of
    // pass 5, F2)
    ms_per_row = std::max(ms_per_row, ms / static_cast<double>(rows));
}

std::vector<BinaryMatrix::SetBitPositions>
LabelOracle::get_rows(const std::vector<Row> &rows) const {
    if (test_read_hook)
        test_read_hook(rows.size());
    Timer timer;
    auto result = matrix_->get_rows(rows);
    counters_.rows_fetched += rows.size();
    counters_.fetch_seconds += timer.elapsed();
    return result;
}

std::vector<MultiIntMatrix::RowTuples>
LabelOracle::get_row_tuples(const std::vector<Row> &rows) const {
    if (!tuples_)
        throw std::logic_error("The annotation has no coordinates");
    if (test_read_hook)
        test_read_hook(rows.size());
    Timer timer;
    auto result = tuples_->get_row_tuples(rows);
    counters_.tuple_rows_fetched += rows.size();
    counters_.fetch_seconds += timer.elapsed();
    return result;
}

bool LabelOracle::get(Row row, Column column) const {
    assert(get_entry_);
    if (test_read_hook)
        test_read_hook(1);
    counters_.direct_reads++;
    return get_entry_->get(row, column);
}

bool LabelOracle::decode_charged() const {
    return rd_ && rd_->supports_budgeted_decode();
}

DecodeStatus LabelOracle::get_rows(const std::vector<Row> &rows, DecodeBudget &budget,
                                   std::vector<BinaryMatrix::SetBitPositions> *out,
                                   std::vector<RowCost> *costs, std::vector<uint64_t> *held) const {
    if (!rd_)
        return DecodeStatus::UNSUPPORTED;
    if (test_read_hook)
        test_read_hook(rows.size());
    Timer timer;
    const DecodeStatus status = rd_->decode_rows(rows, budget, out, costs, held);
    // physical counters (timing): the rows a decode returned, and the time of every decode,
    // a refused one's included — also when the fetch that asked is refused afterwards
    if (status == DecodeStatus::OK)
        counters_.rows_fetched += rows.size();
    counters_.fetch_seconds += timer.elapsed();
    return status;
}

DecodeStatus LabelOracle::get_row_tuples(const std::vector<Row> &rows, DecodeBudget &budget,
                                         std::vector<MultiIntMatrix::RowTuples> *out,
                                         std::vector<RowCost> *costs,
                                         std::vector<uint64_t> *held) const {
    if (!rd_ || !tuples_)
        return DecodeStatus::UNSUPPORTED;
    if (test_read_hook)
        test_read_hook(rows.size());
    Timer timer;
    const DecodeStatus status = rd_->decode_row_tuples(rows, budget, out, costs, held);
    if (status == DecodeStatus::OK)
        counters_.tuple_rows_fetched += rows.size();
    counters_.fetch_seconds += timer.elapsed();
    return status;
}


/********************************* LabelQuery *********************************/

LabelQuery::LabelQuery(const LabelOracle &oracle,
                       std::vector<LabelRef> labels,
                       bool with_coords,
                       LabelOracle::Access access,
                       size_t max_cache_size)
      : oracle_(oracle), labels_(std::move(labels)), with_coords_(with_coords),
        path_(Path::ROWS), max_cache_size_(max_cache_size) {
    if (labels_.empty())
        throw std::invalid_argument("Empty permitted label set");

    bool has_headers = false;
    for (LabelId i = 0; i < labels_.size(); ++i) {
        const LabelRef &ref = labels_[i];
        if (ref.kind == LabelKind::COLUMN) {
            if (!column_labels_.try_emplace(ref.column, i).second)
                throw std::invalid_argument("Duplicate label '" + ref.name + "'");
        } else {
            has_headers = true;
            if (!header_labels_.try_emplace(std::make_pair(ref.column, ref.seq_id), i).second)
                throw std::invalid_argument("Duplicate label '" + ref.name + "'");
            header_columns_.insert(ref.column);
        }
        direct_columns_.push_back(ref.column);
    }
    std::sort(direct_columns_.begin(), direct_columns_.end());
    direct_columns_.erase(std::unique(direct_columns_.begin(), direct_columns_.end()),
                          direct_columns_.end());

    if (has_headers || with_coords_) {
        if (!oracle_.has_coordinates()) {
            throw std::invalid_argument(has_headers
                ? "Sequence header labels require an annotation with coordinates"
                : "Coordinates were requested but the annotation has none");
        }
        if (has_headers && !oracle_.coord_to_header())
            throw std::invalid_argument("Sequence header labels require a CoordToHeader index");
        if (access == LabelOracle::Access::DIRECT)
            throw std::invalid_argument("Direct access is not available for header labels or coordinates");
        path_ = Path::TUPLES;
        return;
    }

    switch (access) {
        case LabelOracle::Access::DIRECT:
            if (!oracle_.supports_direct()) {
                throw std::invalid_argument("Direct access is not supported by this "
                                            "annotation representation");
            }
            path_ = Path::DIRECT;
            break;
        case LabelOracle::Access::ROWS:
            path_ = Path::ROWS;
            break;
        case LabelOracle::Access::AUTO:
            // a few single-cell reads beat reconstructing the whole row
            path_ = oracle_.supports_direct() && direct_columns_.size() <= 16 ? Path::DIRECT
                                                                              : Path::ROWS;
            break;
    }
}

const char* LabelQuery::access_path() const {
    switch (path_) {
        case Path::DIRECT: return "direct";
        case Path::ROWS: return "rows";
        case Path::TUPLES: return "tuples";
    }
    return "unknown";
}

void LabelQuery::hits_from_row(const BinaryMatrix::SetBitPositions &row, NodeHits *hits) const {
    for (Column c : row) {
        auto it = column_labels_.find(c);
        if (it != column_labels_.end())
            hits->push_back(Hit{ it->second, {} });
    }
    std::sort(hits->begin(), hits->end(),
              [](const Hit &a, const Hit &b) { return a.label < b.label; });
}

void LabelQuery::hits_from_tuples(const MultiIntMatrix::RowTuples &row, NodeHits *hits) const {
    // label id -> position in |hits|. A header label is reached once per coordinate,
    // so locating its entry must be O(1): scanning the hits built so far per coordinate
    // made this O(coords x hits) per node, on the coordinate path real indexes use.
    tsl::hopscotch_map<LabelId, size_t> slot;
    for (const auto &[c, coords] : row) {
        if (auto it = column_labels_.find(c); it != column_labels_.end()) {
            Hit hit{ it->second, {} };
            if (with_coords_) {
                hit.coords.assign(coords.begin(), coords.end());
                std::sort(hit.coords.begin(), hit.coords.end());
            }
            slot.emplace(it->second, hits->size());
            hits->push_back(std::move(hit));
        }
        if (!header_columns_.count(c))
            continue;
        for (Coord coord : coords) {
            auto [seq_id, local] = oracle_.map_coord(c, coord);
            auto lt = header_labels_.find(std::make_pair(c, seq_id));
            if (lt == header_labels_.end())
                continue;
            auto [pos, inserted] = slot.try_emplace(lt->second, hits->size());
            if (inserted)
                hits->push_back(Hit{ lt->second, {} });
            if (with_coords_)
                (*hits)[pos->second].coords.push_back(local);
        }
    }
    for (Hit &hit : *hits) {
        std::sort(hit.coords.begin(), hit.coords.end());
    }
    std::sort(hits->begin(), hits->end(),
              [](const Hit &a, const Hit &b) { return a.label < b.label; });
}

void LabelQuery::fetch_uncached(const node_index *keys, size_t n) {
    if (!n)
        return;

    std::vector<Row> rows;
    rows.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        assert(keys[i] != npos);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(keys[i]));
    }

    std::vector<NodeHits> result(n);
    switch (path_) {
        case Path::DIRECT: {
            for (size_t i = 0; i < rows.size(); ++i) {
                for (Column c : direct_columns_) {
                    if (oracle_.get(rows[i], c))
                        result[i].push_back(Hit{ column_labels_.at(c), {} });
                }
                // columns are visited in column order, hits are sorted by label id
                std::sort(result[i].begin(), result[i].end(),
                          [](const Hit &a, const Hit &b) { return a.label < b.label; });
            }
            break;
        }
        case Path::ROWS: {
            auto fetched = oracle_.get_rows(rows);
            for (size_t i = 0; i < rows.size(); ++i) {
                hits_from_row(fetched[i], &result[i]);
            }
            if (max_cache_bytes_ != std::numeric_limits<uint64_t>::max()) {
                // under a byte bound: what the raw rows held beside the hits (observed by
                // the walker as memory_bound_soft), estimated as the cache estimates
                for (const auto &row : fetched) {
                    last_call_bytes_ += sizeof(row) + row.size() * sizeof(row[0]);
                }
            }
            break;
        }
        case Path::TUPLES: {
            auto fetched = oracle_.get_row_tuples(rows);
            for (size_t i = 0; i < rows.size(); ++i) {
                hits_from_tuples(fetched[i], &result[i]);
            }
            if (max_cache_bytes_ != std::numeric_limits<uint64_t>::max()) {
                for (const auto &row : fetched) {
                    last_call_bytes_ += sizeof(row) + row.size() * sizeof(row[0]);
                    for (const auto &entry : row) {
                        last_call_bytes_ += entry.second.size() * sizeof(uint64_t);
                    }
                }
            }
            break;
        }
    }

    // Always insert: a silent no-op here would make the caller's cache_.at() throw.
    // The caller evicts before fetching when the cache would overflow, so the cache can
    // exceed max_cache_size_ only by the working set of a single call.
    for (size_t i = 0; i < n; ++i) {
        // the estimate the byte bound compares: slot, vector, hits and coordinates
        uint64_t bytes = 64 + result[i].size() * sizeof(Hit);
        for (const Hit &h : result[i]) {
            bytes += h.coords.size() * sizeof(Coord);
        }
        cache_bytes_ += bytes;
        cache_[keys[i]] = std::move(result[i]);
    }
}

// |keys| (sorted, distinct) in the order of their first appearance in |order|
static std::vector<node_index> in_first_order(const std::vector<node_index> &keys,
                                              const std::vector<node_index> &order) {
    std::vector<node_index> out;
    out.reserve(keys.size());
    std::vector<bool> taken(keys.size(), false);
    for (node_index key : order) {
        auto it = std::lower_bound(keys.begin(), keys.end(), key);
        if (it != keys.end() && *it == key && !taken[it - keys.begin()]) {
            taken[it - keys.begin()] = true;
            out.push_back(key);
        }
    }
    assert(out.size() == keys.size());
    return out;
}

// The pacing of the unbudgeted reads of LabelQuery and LabelRecorder: |keys| (sorted,
// distinct) decoded by |decode| (their fetch_uncached, which caches what it decodes), every
// piece measured. In one piece, as before pass 5, when the deadline cannot fall into the read
// (DecodePacer::next). Split, its chunks are taken in the order of the keys' first appearance
// in |order| — the caller's keys, in the walk's order, where the nodes of one path are adjacent
// and their rows share the decoding of their row-diff paths, which one call does once — and
// each is sorted as a whole read is: chunks of sorted keys held rows of as many paths as rows,
// and a row_diff lookahead took 0.75 ms a row in chunks against 0.016 ms read whole (review of
// pass 5, F1). Stopped before a chunk: false, ReadPacing::interrupted set and its units the
// work of the keys decoded (|units| of each, read from the cache).
template <class Decode, class Units>
static bool paced_fetch(DecodePacer &pacer, const std::vector<node_index> &keys,
                        const std::vector<node_index> &order, ReadPacing *pacing,
                        const Decode &decode, const Units &units) {
    if (keys.empty())
        return true;
    const bool paced = pacing && pacer.target_ms > 0;
    auto ms_left = [&]() {
        return pacing->ms_left ? pacing->ms_left() : std::numeric_limits<double>::infinity();
    };
    auto stopped = [&](const node_index *decoded, size_t n) {
        if (!pacing->stop || !pacing->stop())
            return false;
        pacing->interrupted = true;
        pacing->units = 0;
        for (size_t i = 0; i < n; ++i) {
            pacing->units += units(decoded[i]);
        }
        return true;
    };
    size_t n = keys.size();
    if (paced) {
        if (stopped(nullptr, 0))
            return false;
        n = pacer.next(keys.size(), ms_left(), 0, 0);
    }
    if (n == keys.size()) {
        Timer timer;
        decode(keys.data(), keys.size());
        pacer.record(keys.size(), timer.elapsed() * 1000);
        return true;
    }
    const std::vector<node_index> ordered = in_first_order(keys, order);
    std::vector<node_index> chunk;
    size_t previous = 0;
    double previous_ms = 0;
    for (size_t begin = 0; begin < ordered.size(); ) {
        if (begin) {
            if (stopped(ordered.data(), begin))
                return false;
            n = pacer.next(ordered.size() - begin, ms_left(), previous, previous_ms);
        }
        chunk.assign(ordered.begin() + begin, ordered.begin() + begin + n);
        std::sort(chunk.begin(), chunk.end());
        Timer timer;
        decode(chunk.data(), n);
        previous_ms = timer.elapsed() * 1000;
        pacer.record(n, previous_ms);
        previous = n;
        begin += n;
    }
    return true;
}

bool LabelQuery::fetch_uncached_paced(const std::vector<node_index> &keys,
                                      const std::vector<node_index> &order, ReadPacing *pacing) {
    return paced_fetch(oracle_.pacer(), keys, order, pacing,
        [this](const node_index *k, size_t n) { fetch_uncached(k, n); },
        [this](node_index key) {
            // the keys decoded so far (in the cache now) are work done: the caller's
            const NodeHits &h = cache_.at(key);
            uint64_t units = 8 + h.size();
            for (const Hit &hit : h) {
                units += hit.coords.size();
            }
            return units;
        });
}

std::vector<LabelQuery::NodeHits> LabelQuery::fetch(const std::vector<node_index> &keys,
                                                    ReadPacing *pacing) {
    last_call_bytes_ = 0;
    // Distinct keys this call must be able to answer. Eviction below may drop entries
    // that were cached on entry, so the whole working set — not just the misses — has
    // to be (re)fetched; otherwise the lookups at the end would throw. The counters are
    // added once the call is answered: a paced call its deadline interrupted changes none
    std::vector<node_index> wanted;
    wanted.reserve(keys.size());
    uint64_t requested = 0, cached = 0;
    for (node_index key : keys) {
        if (key == npos)
            continue;
        requested++;
        if (cache_.count(key))
            cached++;
        wanted.push_back(key);
    }
    std::sort(wanted.begin(), wanted.end());
    wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());

    std::vector<node_index> missing;
    for (node_index key : wanted) {
        if (!cache_.count(key))
            missing.push_back(key);
    }

    if (!missing.empty()) {
        if (cache_.size() + missing.size() > max_cache_size_ || cache_bytes_ > max_cache_bytes_) {
            // a walk moves forward, so evict wholesale rather than tracking recency,
            // then refetch everything this call needs
            clear_cache();
            missing = wanted;
        }
        // decided above for the whole call, so that its chunks decode what one call would
        if (!fetch_uncached_paced(missing, keys, pacing))
            return {};
    }
    oracle_.counters().rows_requested += requested;
    oracle_.counters().cache_hits += cached;

    std::vector<NodeHits> result;
    result.reserve(keys.size());
    for (node_index key : keys) {
        if (key == npos) {
            result.emplace_back();
        } else {
            result.push_back(cache_.at(key));
        }
    }
    return result;
}

void LabelQuery::warm(const std::vector<node_index> &keys, ReadPacing *pacing) {
    last_call_bytes_ = 0;
    std::vector<node_index> missing;
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key))
            missing.push_back(key);
    }
    if (missing.empty())
        return;
    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    if (cache_.size() + missing.size() > max_cache_size_ || cache_bytes_ > max_cache_bytes_)
        clear_cache();
    if (missing.size() >= max_cache_size_)
        return;     // would not fit even into an empty cache: nothing to warm
    fetch_uncached_paced(missing, keys, pacing);
}

const LabelQuery::NodeHits& LabelQuery::fetch(node_index key) {
    if (key == npos)
        return empty_;
    oracle_.counters().rows_requested++;
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        oracle_.counters().cache_hits++;
        return it->second;
    }
    if (cache_.size() >= max_cache_size_ || cache_bytes_ > max_cache_bytes_)
        clear_cache();
    fetch_uncached(&key, 1);
    return cache_.at(key);
}


/********************** LabelQuery: the budget-aware path *********************/

uint64_t LabelQuery::held_bytes(const NodeHits &hits) {
    uint64_t bytes = buffer_bytes(hits.size(), sizeof(Hit));
    for (const Hit &h : hits) {
        bytes += small_vector_bytes(h.coords.size(), sizeof(Coord));
    }
    return bytes;
}

// As hits_from_row, with the hits counted before their vector is reserved at its exact size
bool LabelQuery::hits_budgeted(const BinaryMatrix::SetBitPositions &row, DecodeBudget &budget,
                               NodeHits *hits, uint64_t *peak) {
    size_t m = 0;
    for (Column c : row) {
        m += column_labels_.count(c);
    }
    *peak = buffer_bytes(m, sizeof(Hit));
    if (!budget.charge(*peak))
        return false;
    hits->reserve(m);
    for (Column c : row) {
        auto it = column_labels_.find(c);
        if (it != column_labels_.end())
            hits->push_back(Hit{ it->second, {} });
    }
    std::sort(hits->begin(), hits->end(),
              [](const Hit &a, const Hit &b) { return a.label < b.label; });
    return true;
}

// As hits_from_tuples (the same hits, the same coordinate mappings), in two passes so that
// every buffer is reserved at its exact size and charged before: the first maps the header
// labels' coordinates once (kept beside, then dropped) and counts each label's coordinates,
// the second fills the hits
bool LabelQuery::hits_budgeted(const MultiIntMatrix::RowTuples &row, DecodeBudget &budget,
                               NodeHits *hits, uint64_t *peak) {
    if (label_count_.size() != labels_.size()) {
        // per-label scratch, part of the query's per-label maps
        label_count_.assign(labels_.size(), 0);
        touched_.reserve(labels_.size());
    }
    uint64_t header_coords = 0;
    for (const auto &[c, coords] : row) {
        if (header_columns_.count(c))
            header_coords += coords.size();
    }
    const uint64_t mapped_bytes = buffer_bytes(header_coords, sizeof(std::pair<LabelId, Coord>));
    if (!budget.charge(mapped_bytes))
        return false;
    std::vector<std::pair<LabelId, Coord>> mapped;
    mapped.reserve(header_coords);
    // label_count_[l]: 1 + its coordinates once touched (0: not on this row)
    auto touch = [&](LabelId l, uint64_t coords) {
        if (!label_count_[l]) {
            touched_.push_back(l);
            label_count_[l] = 1;
        }
        label_count_[l] += coords;
    };
    for (const auto &[c, coords] : row) {
        if (auto it = column_labels_.find(c); it != column_labels_.end())
            touch(it->second, with_coords_ ? coords.size() : 0);
        if (!header_columns_.count(c))
            continue;
        for (Coord coord : coords) {
            auto [seq_id, local] = oracle_.map_coord(c, coord);
            auto lt = header_labels_.find(std::make_pair(c, seq_id));
            if (lt == header_labels_.end())
                continue;
            touch(lt->second, with_coords_ ? 1 : 0);
            if (with_coords_)
                mapped.emplace_back(lt->second, local);
        }
    }
    std::sort(touched_.begin(), touched_.end());
    uint64_t hits_bytes = buffer_bytes(touched_.size(), sizeof(Hit));
    for (LabelId l : touched_) {
        hits_bytes += small_vector_bytes(label_count_[l] - 1, sizeof(Coord));
    }
    auto reset = [&]() {
        for (LabelId l : touched_) {
            label_count_[l] = 0;
        }
        touched_.clear();
    };
    if (!budget.charge(hits_bytes)) {
        reset();
        // |mapped| is freed on return
        budget.release(mapped_bytes);
        return false;
    }
    hits->reserve(touched_.size());
    for (size_t i = 0; i < touched_.size(); ++i) {
        const LabelId l = touched_[i];
        hits->push_back(Hit{ l, {} });
        hits->back().coords.reserve(label_count_[l] - 1);
        label_count_[l] = i + 1;     // from here: the label's position in |hits|, plus one
    }
    if (with_coords_) {
        for (const auto &[c, coords] : row) {
            if (auto it = column_labels_.find(c); it != column_labels_.end())
                (*hits)[label_count_[it->second] - 1].coords.assign(coords.begin(), coords.end());
        }
        for (const auto &[l, local] : mapped) {
            (*hits)[label_count_[l] - 1].coords.push_back(local);
        }
        for (Hit &hit : *hits) {
            std::sort(hit.coords.begin(), hit.coords.end());
        }
    }
    reset();
    std::vector<std::pair<LabelId, Coord>>().swap(mapped);
    budget.release(mapped_bytes);
    *peak = mapped_bytes + hits_bytes;
    return true;
}

DecodeStatus LabelQuery::decode_run(const node_index *keys, size_t n, DecodeBudget &budget,
                                    std::vector<NodeHits> *out, size_t at,
                                    std::vector<KeyCost> *costs, size_t *built) {
    *built = 0;
    const uint64_t rows_bytes = buffer_bytes(n, sizeof(Row));
    if (!budget.charge(rows_bytes))
        return DecodeStatus::REFUSED;
    std::vector<Row> rows(n);
    for (size_t i = 0; i < n; ++i) {
        assert(keys[i] != npos);
        rows[i] = AnnotatedDBG::graph_to_anno_index(keys[i]);
    }
    std::vector<RowCost> row_costs;
    std::vector<uint64_t> held;
    std::vector<BinaryMatrix::SetBitPositions> plain;
    std::vector<MultiIntMatrix::RowTuples> tuples;
    const bool tuple_path = path_ == Path::TUPLES;
    assert(path_ != Path::DIRECT);
    const DecodeStatus status = tuple_path
        ? oracle_.get_row_tuples(rows, budget, &tuples, &row_costs, &held)
        : oracle_.get_rows(rows, budget, &plain, &row_costs, &held);
    if (status != DecodeStatus::OK) {
        budget.release(rows_bytes);
        return status;
    }
    const uint64_t containers = tuple_path
        ? annot::matrix::IRowDiff::output_bytes<MultiIntMatrix::RowTuples>(n)
        : annot::matrix::IRowDiff::output_bytes<BinaryMatrix::SetBitPositions>(n);
    // each row's hits replace its raw row; a key's demand: decoding it alone, the row
    // vector of one key, and what building its hits holds beside its raw row
    const uint64_t one = buffer_bytes(1, sizeof(Row));
    DecodeStatus result = DecodeStatus::OK;
    for (size_t i = 0; i < n; ++i) {
        uint64_t peak = 0;
        const bool ok = tuple_path ? hits_budgeted(tuples[i], budget, &(*out)[at + i], &peak)
                                   : hits_budgeted(plain[i], budget, &(*out)[at + i], &peak);
        if (!ok) {
            // the rest of the run is dropped (the caller retries it in a smaller run)
            for (size_t j = i; j < n; ++j) {
                budget.release(held[j]);
            }
            result = DecodeStatus::REFUSED;
            break;
        }
        (*costs)[at + i] = KeyCost{ dependency_units(row_costs[i]),
                                    row_costs[i].demand + one + peak };
        budget.release(held[i]);
        if (tuple_path) {
            MultiIntMatrix::RowTuples().swap(tuples[i]);
        } else {
            BinaryMatrix::SetBitPositions().swap(plain[i]);
        }
        ++*built;
    }
    budget.release(containers + rows_bytes);
    return result;
}

void LabelQuery::cache_budgeted(const node_index *keys, size_t n, const NodeHits *hits,
                                const KeyCost *costs) {
    // What the keys decoded by this call add (|fresh|), and what all its keys would hold in an
    // emptied cache (|all|): a duplicate within the call is counted at each position, an
    // overestimate that needs no call-sized set, whose bytes nothing would charge
    uint64_t fresh = 0, all = 0;
    size_t fresh_count = 0, all_count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (keys[i] == npos)
            continue;
        const uint64_t bytes = kCacheEntryBytes + kCostEntryBytes + held_bytes(hits[i]);
        all += bytes;
        all_count++;
        if (!cache_.count(keys[i])) {
            fresh += bytes;
            fresh_count++;
        }
    }
    if (!fresh_count)
        return;
    // evicted wholesale (a walk moves forward) when they do not fit beside what is cached,
    // keeping this call's keys; not cached at all when they alone exceed the bound: the cache
    // never exceeds it
    if (cache_bytes_ + fresh > max_cache_bytes_ || cache_.size() + fresh_count > max_cache_size_) {
        clear_cache();
        if (all > max_cache_bytes_ || all_count > max_cache_size_)
            return;
    }
    for (size_t i = 0; i < n; ++i) {
        if (keys[i] != npos && !cache_.count(keys[i])) {
            cache_[keys[i]] = hits[i];
            costs_[keys[i]] = costs[i];
            cache_bytes_ += kCacheEntryBytes + kCostEntryBytes + held_bytes(hits[i]);
        }
    }
}

bool LabelQuery::fetch(const node_index *keys, size_t n, DecodeBudget &budget,
                       std::vector<NodeHits> *out, std::vector<KeyCost> *costs,
                       size_t *refused_at, ReadPacing *pacing) {
    // every cached key has its costs (a walk reads by one path; entries the unbudgeted path
    // cached would have none)
    if (costs_.size() != cache_.size())
        clear_cache();
    const uint64_t at_entry = budget.held();
    const size_t base = out->size();
    assert(costs->size() == base);
    assert(out->capacity() >= base + n && costs->capacity() >= base + n);
    // what the keys admitted so far hold (their hits); the cache is not changed before the
    // end, so the keys it does not hold then are the ones this call decoded
    uint64_t committed = 0;
    size_t requested = 0, hits_from_cache = 0;
    auto left = [&]() { return budget.max_bytes() - at_entry - committed; };
    auto refuse = [&](size_t pos, FetchRefusal::Cause cause, uint64_t demand, uint64_t need) {
        refusal_ = FetchRefusal();
        refusal_.cause = cause;
        refusal_.position = pos;
        refusal_.left = left();
        refusal_.held = committed;
        refusal_.demand = demand;
        refusal_.need = need;
        out->resize(base);
        costs->resize(base);
        budget.restore(at_entry);
        *refused_at = pos;
        return false;
    };
    // the least a refused read alone was seen to need: what it held when its charge did not
    // fit, and that charge
    auto seen_need = [&]() {
        const uint64_t before = at_entry + committed;
        return budget.refused_need() > before ? budget.refused_need() - before : 0;
    };
    // A paced fetch's deadline before a run: restored as a refusal (nothing appended, the
    // budget as on entry), with the work of the keys its runs decoded (those the cache, which
    // does not change during the call, does not hold) for the caller to charge
    DecodePacer &pacer = oracle_.pacer();
    const bool paced = pacing && pacer.target_ms > 0;
    size_t previous = 0;
    double previous_ms = 0;
    auto interrupt = [&](size_t pos) {
        pacing->interrupted = true;
        pacing->units = 0;
        for (size_t i = 0; i < pos; ++i) {
            if (keys[i] == npos || cache_.count(keys[i]))
                continue;
            const NodeHits &h = (*out)[base + i];
            pacing->units += 8 + h.size() + (*costs)[base + i].dependency_units;
            for (const Hit &hit : h) {
                pacing->units += hit.coords.size();
            }
        }
        refuse(pos, FetchRefusal::INTERRUPTED, 0, 0);
        return false;
    };
    size_t run_limit = kMaxDecodeRun;
    for (size_t pos = 0; pos < n; ) {
        const node_index key = keys[pos];
        if (key == npos) {
            out->emplace_back();
            costs->emplace_back();
            ++pos;
            continue;
        }
        auto it = cache_.find(key);
        if (it != cache_.end()) {
            // admitted against its demand, as if decoded now (its costs are the key's)
            const KeyCost &cost = costs_.at(key);
            if (cost.demand > left())
                return refuse(pos, FetchRefusal::DEMAND, cost.demand, 0);
            // the demand covers the copy, so only the test hook can refuse it
            const uint64_t bytes = held_bytes(it->second);
            if (!budget.charge(bytes))
                return refuse(pos, FetchRefusal::DECODE, 0, seen_need());
            committed += bytes;
            out->push_back(it->second);
            costs->push_back(cost);
            requested++;
            hits_from_cache++;
            ++pos;
            continue;
        }
        // a run of consecutive misses, decoded together (at most a paced chunk: admission is
        // per key, so how the misses are cut into runs changes nothing that is returned)
        size_t limit = run_limit;
        if (paced) {
            if (pacing->stop && pacing->stop())
                return interrupt(pos);
            limit = std::min(limit, pacer.next(n - pos, pacing->ms_left
                                                            ? pacing->ms_left()
                                                            : std::numeric_limits<double>::infinity(),
                                               previous, previous_ms));
        }
        size_t end = pos;
        while (end < n && end - pos < limit && keys[end] != npos && !cache_.count(keys[end])) {
            ++end;
        }
        const size_t len = end - pos;
        out->resize(base + end);
        costs->resize(base + end);
        size_t built = 0;
        Timer timer;
        const DecodeStatus status = decode_run(keys + pos, len, budget, out, base + pos, costs,
                                               &built);
        previous_ms = timer.elapsed() * 1000;
        pacer.record(len, previous_ms);
        previous = len;
        // the keys built, each admitted against what was left at its position
        for (size_t j = 0; j < built; ++j) {
            const uint64_t demand = (*costs)[base + pos + j].demand;
            if (demand > left())
                return refuse(pos + j, FetchRefusal::DEMAND, demand, 0);
            committed += held_bytes((*out)[base + pos + j]);
        }
        assert(at_entry + committed == budget.held());
        requested += built;
        if (status == DecodeStatus::OK) {
            pos = end;
            continue;
        }
        assert(status == DecodeStatus::REFUSED);
        // the run did not fit as one: a single key that does not fit alone is the stop;
        // otherwise go on from the first key not built, in smaller runs
        if (!built && len == 1)
            return refuse(pos, FetchRefusal::DECODE, 0, seen_need());
        pos += built;
        out->resize(base + pos);
        costs->resize(base + pos);
        run_limit = std::max<size_t>(1, len / 2);
    }
    oracle_.counters().rows_requested += requested;
    oracle_.counters().cache_hits += hits_from_cache;
    cache_budgeted(keys, n, out->data() + base, costs->data() + base);
    return true;
}

void LabelQuery::warm(const std::vector<node_index> &keys, DecodeBudget &budget,
                      ReadPacing *pacing) {
    if (costs_.size() != cache_.size())
        clear_cache();
    // A cache of capacity zero keeps nothing, so there is nothing to warm; the runs below
    // are at most that capacity long and would never advance (review of stage 3, F1: a
    // warm with max_cache_size 0 did not return)
    if (!max_cache_size_)
        return;
    // the misses, charged like the rest of the lookahead's read
    const uint64_t at_entry = budget.held();
    if (!budget.charge(buffer_bytes(keys.size(), sizeof(node_index))))
        return;
    std::vector<node_index> missing;
    missing.reserve(keys.size());
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key))
            missing.push_back(key);
    }
    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    DecodePacer &pacer = oracle_.pacer();
    const bool paced = pacing && pacer.target_ms > 0;
    size_t previous = 0;
    double previous_ms = 0;
    const size_t run = std::min(kMaxDecodeRun, max_cache_size_);
    for (size_t begin = 0; begin < missing.size(); begin += run) {
        const size_t len = std::min(run, missing.size() - begin);
        const uint64_t at_run = budget.held();
        const uint64_t containers = buffer_bytes(len, sizeof(NodeHits))
                                  + buffer_bytes(len, sizeof(KeyCost));
        if (!budget.charge(containers))
            break;
        std::vector<NodeHits> hits(len);
        std::vector<KeyCost> costs(len);
        // the run decoded in paced pieces, cached as one run (the cache's decisions are the
        // run's); a deadline before a piece drops the run and ends the warming
        bool ok = true;
        for (size_t at = 0; ok && at < len; ) {
            size_t piece = len - at;
            if (paced) {
                if (pacing->stop && pacing->stop()) {
                    pacing->interrupted = true;
                    budget.restore(at_entry);
                    return;
                }
                piece = pacer.next(piece, pacing->ms_left ? pacing->ms_left()
                                                          : std::numeric_limits<double>::infinity(),
                                   previous, previous_ms);
            }
            size_t built = 0;
            Timer timer;
            ok = decode_run(missing.data() + begin + at, piece, budget, &hits, at, &costs,
                            &built) == DecodeStatus::OK;
            previous_ms = timer.elapsed() * 1000;
            pacer.record(piece, previous_ms);
            previous = piece;
            at += piece;
        }
        if (!ok)
            break;      // the lookahead gives up within what is left: nothing depends on it
        cache_budgeted(missing.data() + begin, len, hits.data(), costs.data());
        budget.restore(at_run);
    }
    budget.restore(at_entry);
}


/******************************** LabelRecorder *******************************/

LabelRecorder::LabelRecorder(const LabelOracle &oracle,
                             LabelKind kind,
                             size_t max_labels_per_node,
                             size_t max_cache_size,
                             size_t max_cache_keys)
      : oracle_(oracle), kind_(kind), cap_(max_labels_per_node),
        max_cache_size_(max_cache_size), max_cache_keys_(max_cache_keys) {
    if (!cap_)
        throw std::invalid_argument("max_labels_per_node must be positive");
    if (kind_ == LabelKind::HEADER) {
        if (!oracle_.has_coordinates()) {
            throw std::invalid_argument("Recording sequence header labels requires an "
                                        "annotation with coordinates");
        }
        if (!oracle_.coord_to_header()) {
            throw std::invalid_argument("Recording sequence header labels requires a "
                                        "CoordToHeader index; set seed_label_kind to \"column\"");
        }
    }
}

const char* LabelRecorder::access_path() const {
    return kind_ == LabelKind::HEADER ? "tuples" : "rows";
}

std::optional<LabelId> LabelRecorder::named(const Key &key) const {
    if (kind_ == LabelKind::COLUMN) {
        auto it = column_ids_.find(key.first);
        if (it != column_ids_.end())
            return it->second;
    } else {
        auto it = header_ids_.find(key);
        if (it != header_ids_.end())
            return it->second;
    }
    return std::nullopt;
}

LabelId LabelRecorder::id_of(const Key &key) {
    if (auto id = named(key))
        return *id;
    LabelRef ref;
    ref.kind = kind_;
    ref.column = key.first;
    ref.seq_id = kind_ == LabelKind::HEADER ? key.second : 0;
    ref.name = kind_ == LabelKind::COLUMN ? oracle_.column_name(key.first)
                                          : oracle_.header_name(key.first, key.second);
    const LabelId id = dict_.size();
    dict_.push_back(std::move(ref));
    if (kind_ == LabelKind::COLUMN) {
        column_ids_.emplace(key.first, id);
    } else {
        header_ids_.emplace(key, id);
    }
    return id;
}

void LabelRecorder::fetch_uncached(const node_index *keys, size_t n) {
    if (!n)
        return;
    std::vector<Row> rows;
    rows.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        assert(keys[i] != npos);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(keys[i]));
    }
    std::vector<RawRow> result(n);
    const bool bounded = max_cache_bytes_ != std::numeric_limits<uint64_t>::max();
    if (kind_ == LabelKind::COLUMN) {
        auto fetched = oracle_.get_rows(rows);
        for (size_t i = 0; i < rows.size(); ++i) {
            // SetBitPositions are ascending columns
            RawRow &r = result[i];
            r.total = fetched[i].size();
            for (size_t j = 0; j < fetched[i].size() && j < cap_; ++j) {
                r.kept.emplace_back(fetched[i][j], 0);
            }
            // under a byte bound: what the raw row held (observed as memory_bound_soft)
            if (bounded)
                last_call_bytes_ += sizeof(fetched[i]) + fetched[i].size() * sizeof(Column);
        }
    } else {
        auto fetched = oracle_.get_row_tuples(rows);
        if (bounded) {
            for (const auto &row : fetched) {
                last_call_bytes_ += sizeof(row) + row.size() * sizeof(row[0]);
                for (const auto &entry : row) {
                    last_call_bytes_ += entry.second.size() * sizeof(uint64_t);
                }
            }
        }
        std::vector<Key> all;
        std::vector<Coord> coords;
        for (size_t i = 0; i < rows.size(); ++i) {
            all.clear();
            for (const auto &[c, tuple] : fetched[i]) {
                // The true count needs to know which SEQUENCES of the column this k-mer
                // belongs to, and only the coordinate mapping tells. But the sequences
                // of a column occupy contiguous coordinate ranges, so one mapping per
                // sequence is enough: map the first coordinate of a sequence, then skip
                // to the first coordinate past that sequence's range. A k-mer repeated
                // a thousand times in one record costs one rank/select, not a thousand.
                coords.assign(tuple.begin(), tuple.end());
                std::sort(coords.begin(), coords.end());
                for (auto it = coords.begin(); it != coords.end(); ) {
                    const auto [seq_id, local] = oracle_.map_coord(c, *it);
                    all.emplace_back(c, seq_id);
                    const Coord end = *it - local + oracle_.num_kmers_in_sequence(c, seq_id);
                    it = std::lower_bound(it, coords.end(), end);
                }
            }
            std::sort(all.begin(), all.end());
            all.erase(std::unique(all.begin(), all.end()), all.end());
            RawRow &r = result[i];
            r.total = all.size();
            r.kept.assign(all.begin(), all.begin() + std::min(all.size(), cap_));
        }
    }
    for (size_t i = 0; i < n; ++i) {
        // a row costs at least one slot even when nothing is on it
        cached_keys_ += std::max<size_t>(1, result[i].kept.size());
        cache_bytes_ += 64 + result[i].kept.size() * sizeof(Key);
        cache_[keys[i]] = std::move(result[i]);
    }
}

bool LabelRecorder::fetch_uncached_paced(const std::vector<node_index> &keys,
                                         const std::vector<node_index> &order,
                                         ReadPacing *pacing) {
    return paced_fetch(oracle_.pacer(), keys, order, pacing,
        [this](const node_index *k, size_t n) { fetch_uncached(k, n); },
        // a recorded row is its whole width: the true count, not the capped list
        [this](node_index key) { return 8 + static_cast<uint64_t>(cache_.at(key).total); });
}

std::vector<LabelRecorder::NodeLabels>
LabelRecorder::fetch(const std::vector<node_index> &keys, ReadPacing *pacing) {
    last_call_bytes_ = 0;
    // as LabelQuery::fetch: the counters are added, and the labels named, once the call is
    // answered, so that an interrupted call changes neither
    std::vector<node_index> wanted;
    wanted.reserve(keys.size());
    uint64_t requested = 0, cached = 0;
    for (node_index key : keys) {
        if (key == npos)
            continue;
        requested++;
        if (cache_.count(key))
            cached++;
        wanted.push_back(key);
    }
    std::sort(wanted.begin(), wanted.end());
    wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());
    std::vector<node_index> missing;
    for (node_index key : wanted) {
        if (!cache_.count(key))
            missing.push_back(key);
    }
    if (!missing.empty()) {
        // the cache is bounded in rows AND in kept keys (its memory), since a row holds
        // up to |cap_| of them; either bound exceeded evicts wholesale, and the whole
        // working set is then refetched so that the lookups below cannot throw
        if (cache_.size() + missing.size() > max_cache_size_ || cached_keys_ > max_cache_keys_
                || cache_bytes_ > max_cache_bytes_) {
            clear_cache();
            missing = wanted;
        }
        if (!fetch_uncached_paced(missing, keys, pacing))
            return {};
    }
    oracle_.counters().rows_requested += requested;
    oracle_.counters().cache_hits += cached;
    std::vector<NodeLabels> result;
    result.reserve(keys.size());
    for (node_index key : keys) {
        NodeLabels nl;
        if (key != npos) {
            const RawRow &raw = cache_.at(key);
            nl.total = raw.total;
            nl.labels.reserve(raw.kept.size());
            // ids are assigned in consumption order; the list is reported by id
            for (const Key &k : raw.kept) {
                nl.labels.push_back(id_of(k));
            }
            std::sort(nl.labels.begin(), nl.labels.end());
        }
        result.push_back(std::move(nl));
    }
    return result;
}

void LabelRecorder::warm(const std::vector<node_index> &keys, ReadPacing *pacing) {
    last_call_bytes_ = 0;
    std::vector<node_index> missing;
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key))
            missing.push_back(key);
    }
    if (missing.empty())
        return;
    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    if (cache_.size() + missing.size() > max_cache_size_ || cached_keys_ > max_cache_keys_
            || cache_bytes_ > max_cache_bytes_) {
        clear_cache();
    }
    if (missing.size() >= max_cache_size_)
        return;
    fetch_uncached_paced(missing, keys, pacing);
}


/******************** LabelRecorder: the budget-aware path ********************/

namespace {

// The labels a budget-aware LabelRecorder::fetch names provisionally — given to the dictionary
// only once the whole call is admitted — with their ids, in the order the call meets them: an
// open-addressing table at load at most 1/2 and a list, both grown by fixed policies (16
// entries, then doubling) and nothing else, so that what m labels make them hold at once,
// the transients of a rehash and of a copy included, is the function pending_bytes(m) of
// LabelRecorder. (A map per column, as before, held a whole bucket array per new label,
// uncharged: review of stage 3, F1.)
class PendingLabels {
  public:
    using Key = std::pair<Column, uint64_t>;

    // a seq_id no label has: the empty slot
    static constexpr uint64_t kEmpty = std::numeric_limits<uint64_t>::max();

    std::optional<LabelId> find(const Key &key) const {
        if (slots_.empty())
            return std::nullopt;
        for (size_t i = slot_of(key); ; i = (i + 1) & (slots_.size() - 1)) {
            if (slots_[i].key == key)
                return slots_[i].id;
            if (slots_[i].key.second == kEmpty)
                return std::nullopt;
        }
    }
    // |key| must not be pending
    void add(const Key &key, LabelId id) {
        assert(key.second != kEmpty && !find(key));
        if (order_.size() == order_.capacity())
            order_.reserve(std::max<size_t>(16, 2 * order_.capacity()));
        order_.push_back(key);
        if (2 * order_.size() > slots_.size()) {
            std::vector<Slot> slots(std::max<size_t>(16, 2 * slots_.size()));
            slots.swap(slots_);
            for (const Slot &slot : slots) {
                if (slot.key.second != kEmpty)
                    place(slot);
            }
        }
        place(Slot{ key, id });
    }
    const std::vector<Key>& order() const { return order_; }

    // the most |m| labels make the table and the list hold at once: their last buffers, and
    // during a growth the one before
    static uint64_t bytes_for(uint64_t m) {
        uint64_t slots = 0, prev_slots = 0, list = 0, prev_list = 0;
        for (uint64_t size = 1; size <= m; ++size) {
            if (size > list) {
                prev_list = list;
                list = std::max<uint64_t>(16, 2 * list);
            }
            if (2 * size > slots) {
                prev_slots = slots;
                slots = std::max<uint64_t>(16, 2 * slots);
            }
        }
        return buffer_bytes(slots, sizeof(Slot)) + buffer_bytes(prev_slots, sizeof(Slot))
                + buffer_bytes(list, sizeof(Key)) + buffer_bytes(prev_list, sizeof(Key));
    }

  private:
    struct Slot {
        Key key { 0, kEmpty };
        LabelId id = 0;
    };

    size_t slot_of(const Key &key) const {
        return LabelKeyHash()(key) & (slots_.size() - 1);
    }
    void place(const Slot &slot) {
        size_t i = slot_of(slot.key);
        while (slots_[i].key.second != kEmpty) {
            i = (i + 1) & (slots_.size() - 1);
        }
        slots_[i] = slot;
    }

    std::vector<Slot> slots_;
    std::vector<Key> order_;
};

} // namespace

uint64_t LabelRecorder::pending_bytes(uint64_t m) {
    return PendingLabels::bytes_for(m);
}

uint64_t LabelRecorder::held_bytes(const NodeLabels &labels) {
    return buffer_bytes(labels.labels.size(), sizeof(LabelId));
}

uint64_t LabelRecorder::raw_bytes(const RawRow &raw) {
    return buffer_bytes(raw.kept.size(), sizeof(Key));
}

// a run's raw rows, their held bytes and their costs while its keys are listed
static uint64_t run_bytes_of(size_t n, size_t raw_row_size) {
    return buffer_bytes(n, raw_row_size) + buffer_bytes(n, sizeof(uint64_t))
            + buffer_bytes(n, sizeof(KeyCost));
}

// As fetch_uncached's COLUMN branch: the first |cap_| columns, reserved at their number
bool LabelRecorder::raw_budgeted(const BinaryMatrix::SetBitPositions &row, DecodeBudget &budget,
                                 RawRow *raw, uint64_t *peak) {
    const size_t kept = std::min<size_t>(row.size(), cap_);
    *peak = buffer_bytes(kept, sizeof(Key));
    if (!budget.charge(*peak))
        return false;
    raw->total = row.size();
    raw->kept.reserve(kept);
    for (size_t j = 0; j < kept; ++j) {
        raw->kept.emplace_back(row[j], 0);
    }
    return true;
}

// As fetch_uncached's HEADER branch (one coordinate mapping per sequence a column's
// coordinates enter), with its two working vectors reserved at their largest sizes (every
// coordinate its own sequence; the widest tuple) and charged before
bool LabelRecorder::raw_budgeted(const MultiIntMatrix::RowTuples &row, DecodeBudget &budget,
                                 RawRow *raw, uint64_t *peak) {
    uint64_t coords_total = 0, widest = 0;
    for (const auto &entry : row) {
        coords_total += entry.second.size();
        widest = std::max<uint64_t>(widest, entry.second.size());
    }
    const uint64_t work = buffer_bytes(coords_total, sizeof(Key)) + buffer_bytes(widest, sizeof(Coord));
    if (!budget.charge(work))
        return false;
    std::vector<Key> all;
    std::vector<Coord> coords;
    all.reserve(coords_total);
    coords.reserve(widest);
    for (const auto &[c, tuple] : row) {
        coords.assign(tuple.begin(), tuple.end());
        std::sort(coords.begin(), coords.end());
        for (auto it = coords.begin(); it != coords.end(); ) {
            const auto [seq_id, local] = oracle_.map_coord(c, *it);
            all.emplace_back(c, seq_id);
            const Coord end = *it - local + oracle_.num_kmers_in_sequence(c, seq_id);
            it = std::lower_bound(it, coords.end(), end);
        }
    }
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    const size_t kept = std::min(all.size(), cap_);
    const uint64_t kept_bytes = buffer_bytes(kept, sizeof(Key));
    if (!budget.charge(kept_bytes)) {
        budget.release(work);
        return false;
    }
    raw->total = all.size();
    raw->kept.assign(all.begin(), all.begin() + kept);
    std::vector<Key>().swap(all);
    std::vector<Coord>().swap(coords);
    budget.release(work);
    *peak = work + kept_bytes;
    return true;
}

DecodeStatus LabelRecorder::decode_run(const node_index *keys, size_t n, DecodeBudget &budget,
                                       std::vector<RawRow> *rows_out,
                                       std::vector<uint64_t> *rows_held,
                                       std::vector<KeyCost> *costs, size_t *built, size_t at) {
    *built = 0;
    const uint64_t rows_bytes = buffer_bytes(n, sizeof(Row));
    if (!budget.charge(rows_bytes))
        return DecodeStatus::REFUSED;
    std::vector<Row> rows(n);
    for (size_t i = 0; i < n; ++i) {
        assert(keys[i] != npos);
        rows[i] = AnnotatedDBG::graph_to_anno_index(keys[i]);
    }
    std::vector<RowCost> row_costs;
    std::vector<uint64_t> held;
    std::vector<BinaryMatrix::SetBitPositions> plain;
    std::vector<MultiIntMatrix::RowTuples> tuples;
    const bool tuple_path = kind_ == LabelKind::HEADER;
    const DecodeStatus status = tuple_path
        ? oracle_.get_row_tuples(rows, budget, &tuples, &row_costs, &held)
        : oracle_.get_rows(rows, budget, &plain, &row_costs, &held);
    if (status != DecodeStatus::OK) {
        budget.release(rows_bytes);
        return status;
    }
    const uint64_t containers = tuple_path
        ? annot::matrix::IRowDiff::output_bytes<MultiIntMatrix::RowTuples>(n)
        : annot::matrix::IRowDiff::output_bytes<BinaryMatrix::SetBitPositions>(n);
    // the row vector and the run's vectors of one key (fetch())
    const uint64_t one = buffer_bytes(1, sizeof(Row)) + run_bytes_of(1, sizeof(RawRow));
    DecodeStatus result = DecodeStatus::OK;
    for (size_t i = 0; i < n; ++i) {
        uint64_t peak = 0;
        const bool ok = tuple_path ? raw_budgeted(tuples[i], budget, &(*rows_out)[at + i], &peak)
                                   : raw_budgeted(plain[i], budget, &(*rows_out)[at + i], &peak);
        if (!ok) {
            for (size_t j = i; j < n; ++j) {
                budget.release(held[j]);
            }
            result = DecodeStatus::REFUSED;
            break;
        }
        // a key's demand: decoding it alone, the row vector of one key, building its raw
        // row, and the label list it is returned as
        (*costs)[at + i] = KeyCost{ dependency_units(row_costs[i]),
                                    row_costs[i].demand + one + peak
                                        + buffer_bytes((*rows_out)[at + i].kept.size(),
                                                       sizeof(LabelId)) };
        (*rows_held)[at + i] = raw_bytes((*rows_out)[at + i]);
        budget.release(held[i]);
        if (tuple_path) {
            MultiIntMatrix::RowTuples().swap(tuples[i]);
        } else {
            BinaryMatrix::SetBitPositions().swap(plain[i]);
        }
        ++*built;
    }
    budget.release(containers + rows_bytes);
    return result;
}

void LabelRecorder::cache_raw(node_index key, RawRow &&raw, const KeyCost &cost) {
    cached_keys_ += std::max<size_t>(1, raw.kept.size());
    cache_bytes_ += kCacheEntryBytes + kCostEntryBytes + raw_bytes(raw);
    cache_[key] = std::move(raw);
    costs_[key] = cost;
}

bool LabelRecorder::fetch(const node_index *keys, size_t n, DecodeBudget &budget,
                          std::vector<NodeLabels> *out, std::vector<KeyCost> *costs,
                          size_t *refused_at,
                          const std::function<uint64_t(std::string_view name)> &name_bytes,
                          ReadPacing *pacing) {
    if (costs_.size() != cache_.size())
        clear_cache();
    const uint64_t at_entry = budget.held();
    const size_t base = out->size();
    assert(costs->size() == base);
    assert(out->capacity() >= base + n && costs->capacity() >= base + n);
    // The labels this call names, in the order fetch() would name them, with their ids
    // continuing the dictionary: given (id_of) only once the whole call is admitted. Each new
    // label is admitted with its priced name and kNamingBytes, which bound what |pending|
    // holds (pending_bytes), so that its charge does not depend on the labels before it.
    PendingLabels pending;
    auto known = [&](const Key &k) -> std::optional<LabelId> {
        if (auto id = named(k))
            return id;
        return pending.find(k);
    };
    // what naming the new labels of |raw| costs: their names as the caller's account will
    // charge them, read in place (a copy of a name would be held uncharged), and the naming
    auto names_of = [&](const RawRow &raw, uint64_t *labels) {
        uint64_t bytes = 0;
        *labels = 0;
        // the kept keys of a row are distinct
        for (const Key &k : raw.kept) {
            if (known(k))
                continue;
            const std::string &name = kind_ == LabelKind::COLUMN
                ? oracle_.column_name(k.first) : oracle_.header_name(k.first, k.second);
            bytes += name_bytes(name) + kNamingBytes;
            ++*labels;
        }
        return bytes;
    };
    // the list of |raw|, naming its new labels (provisionally)
    auto list_of = [&](const RawRow &raw, NodeLabels *nl) {
        nl->total = raw.total;
        nl->labels.reserve(raw.kept.size());
        for (const Key &k : raw.kept) {
            auto id = known(k);
            if (!id) {
                id = static_cast<LabelId>(dict_.size() + pending.order().size());
                pending.add(k, *id);
            }
            nl->labels.push_back(*id);
        }
        std::sort(nl->labels.begin(), nl->labels.end());
    };
    // what the keys admitted so far hold: their lists, the names they gave and their naming
    uint64_t committed = 0, names_given = 0, naming = 0;
    auto left = [&]() { return budget.max_bytes() - at_entry - committed; };
    auto refuse = [&](size_t pos, FetchRefusal::Cause cause, uint64_t demand, uint64_t need,
                      uint64_t labels, uint64_t names) {
        refusal_ = FetchRefusal();
        refusal_.cause = cause;
        refusal_.position = pos;
        refusal_.left = left();
        refusal_.held = committed;
        refusal_.demand = demand;
        refusal_.need = need;
        refusal_.labels = labels;
        refusal_.names_bytes = names;
        out->resize(base);
        costs->resize(base);
        budget.restore(at_entry);
        *refused_at = pos;
        return false;
    };
    auto seen_need = [&]() {
        const uint64_t before = at_entry + committed;
        return budget.refused_need() > before ? budget.refused_need() - before : 0;
    };
    // a key read (cached or decoded now) is admitted against its demand and its names, and a
    // refusal says which of them did not fit
    auto admit = [&](size_t pos, uint64_t demand, uint64_t labels, uint64_t names) {
        if (demand > left())
            return refuse(pos, FetchRefusal::DEMAND, demand, 0, labels, names);
        if (demand + names > left())
            return refuse(pos, FetchRefusal::NAMES, demand, 0, labels, names);
        return true;
    };
    // as LabelQuery's: a paced fetch's deadline before a run is restored as a refusal (and
    // names nothing), with the work of the keys its runs decoded for the caller to charge
    DecodePacer &pacer = oracle_.pacer();
    const bool paced = pacing && pacer.target_ms > 0;
    size_t previous = 0;
    double previous_ms = 0;
    auto interrupt = [&](size_t pos) {
        pacing->interrupted = true;
        pacing->units = 0;
        for (size_t i = 0; i < pos; ++i) {
            if (keys[i] != npos && !cache_.count(keys[i]))
                pacing->units += 8 + (*out)[base + i].total + (*costs)[base + i].dependency_units;
        }
        return refuse(pos, FetchRefusal::INTERRUPTED, 0, 0, 0, 0);
    };
    size_t requested = 0, hits_from_cache = 0;
    size_t run_limit = kMaxDecodeRun;
    for (size_t pos = 0; pos < n; ) {
        const node_index key = keys[pos];
        if (key == npos) {
            out->emplace_back();
            costs->emplace_back();
            ++pos;
            continue;
        }
        auto it = cache_.find(key);
        if (it != cache_.end()) {
            const KeyCost &cost = costs_.at(key);
            uint64_t labels = 0;
            const uint64_t names = names_of(it->second, &labels);
            if (!admit(pos, cost.demand, labels, names))
                return false;
            NodeLabels nl;
            // the demand covers the list, so only the test hook can refuse it
            const uint64_t list = buffer_bytes(it->second.kept.size(), sizeof(LabelId));
            if (!budget.charge(list + names))
                return refuse(pos, FetchRefusal::DECODE, 0, seen_need(), 0, 0);
            committed += list + names;
            naming += labels * kNamingBytes;
            names_given += names - labels * kNamingBytes;
            list_of(it->second, &nl);
            out->push_back(std::move(nl));
            costs->push_back(cost);
            requested++;
            hits_from_cache++;
            ++pos;
            continue;
        }
        // a run of consecutive misses, decoded together into raw rows (the cache's form),
        // then listed key by key; at most a paced chunk long
        size_t limit = run_limit;
        if (paced) {
            if (pacing->stop && pacing->stop())
                return interrupt(pos);
            limit = std::min(limit, pacer.next(n - pos, pacing->ms_left
                                                            ? pacing->ms_left()
                                                            : std::numeric_limits<double>::infinity(),
                                               previous, previous_ms));
        }
        size_t end = pos;
        while (end < n && end - pos < limit && keys[end] != npos && !cache_.count(keys[end])) {
            ++end;
        }
        const size_t len = end - pos;
        const uint64_t work = run_bytes_of(len, sizeof(RawRow));
        if (!budget.charge(work)) {
            if (len == 1)
                return refuse(pos, FetchRefusal::DECODE, 0, seen_need(), 0, 0);
            run_limit = std::max<size_t>(1, len / 2);
            continue;
        }
        std::vector<RawRow> raws(len);
        std::vector<uint64_t> raws_held(len);
        std::vector<KeyCost> run_costs(len);
        size_t built = 0;
        Timer timer;
        DecodeStatus status = decode_run(keys + pos, len, budget, &raws, &raws_held, &run_costs,
                                         &built);
        previous_ms = timer.elapsed() * 1000;
        pacer.record(len, previous_ms);
        previous = len;
        if (!built && len == 1) {
            assert(status == DecodeStatus::REFUSED);
            // seen before the run's own buffers are freed: they are part of the read
            return refuse(pos, FetchRefusal::DECODE, 0, seen_need(), 0, 0);
        }
        size_t listed = 0;
        for (; listed < built; ++listed) {
            const RawRow &raw = raws[listed];
            uint64_t labels = 0;
            const uint64_t names = names_of(raw, &labels);
            // admitted against what was left at its position, as a cached key is: the run's
            // buffers are held beside it (charged), which its demand counts for one key
            if (!admit(pos + listed, run_costs[listed].demand, labels, names))
                return false;
            const uint64_t list = buffer_bytes(raw.kept.size(), sizeof(LabelId));
            if (!budget.charge(list + names)) {
                // the rest of the run held too much beside it: go on in smaller runs
                status = DecodeStatus::REFUSED;
                break;
            }
            NodeLabels nl;
            list_of(raw, &nl);
            out->push_back(std::move(nl));
            costs->push_back(run_costs[listed]);
            committed += list + names;
            naming += labels * kNamingBytes;
            names_given += names - labels * kNamingBytes;
            budget.release(raws_held[listed]);
            RawRow().kept.swap(raws[listed].kept);
            requested++;
        }
        for (size_t j = listed; j < built; ++j) {
            budget.release(raws_held[j]);
        }
        raws = std::vector<RawRow>();
        raws_held = std::vector<uint64_t>();
        run_costs = std::vector<KeyCost>();
        budget.release(work);
        if (status == DecodeStatus::OK) {
            pos = end;
            continue;
        }
        if (!listed && len == 1)
            return refuse(pos, FetchRefusal::DECODE, 0, seen_need(), 0, 0);
        pos += listed;
        run_limit = std::max<size_t>(1, len / 2);
    }
    // admitted whole: the labels are named, in the order they were met
    assert(pending.order().size() * kNamingBytes == naming);
    assert(PendingLabels::bytes_for(pending.order().size()) <= naming || !naming);
    for (const Key &k : pending.order()) {
        const LabelId id = id_of(k);
        assert(id + 1 == dict_.size());
        (void)id;
    }
    oracle_.counters().rows_requested += requested;
    oracle_.counters().cache_hits += hits_from_cache;
    last_names_bytes_ = names_given;
    last_naming_bytes_ = naming;
    // The keys this call decoded (those the cache does not hold: it has not changed during the
    // call), cached with their costs within the bounds, their raw rows rebuilt from the lists
    // (a list is the ids of the raw row's kept labels). A duplicate within the call is counted
    // at each position: an overestimate that needs no call-sized set, whose bytes nothing
    // would charge. They do not fit beside what is cached: evicted wholesale, keeping this
    // call's keys; they alone exceed the bounds: not cached.
    uint64_t fresh = 0, all = 0, fresh_kept = 0, all_kept = 0;
    size_t fresh_count = 0, all_count = 0;
    for (size_t p = 0; p < n; ++p) {
        if (keys[p] == npos)
            continue;
        const NodeLabels &nl = (*out)[base + p];
        const uint64_t bytes = kCacheEntryBytes + kCostEntryBytes
                             + buffer_bytes(nl.labels.size(), sizeof(Key));
        const uint64_t kept = std::max<size_t>(1, nl.labels.size());
        all += bytes;
        all_kept += kept;
        all_count++;
        if (!cache_.count(keys[p])) {
            fresh += bytes;
            fresh_kept += kept;
            fresh_count++;
        }
    }
    if (!fresh_count)
        return true;
    if (cache_bytes_ + fresh > max_cache_bytes_ || cache_.size() + fresh_count > max_cache_size_
            || cached_keys_ + fresh_kept > max_cache_keys_) {
        clear_cache();
        if (all > max_cache_bytes_ || all_count > max_cache_size_ || all_kept > max_cache_keys_)
            return true;
    }
    for (size_t p = 0; p < n; ++p) {
        if (keys[p] == npos || cache_.count(keys[p]))
            continue;
        const NodeLabels &nl = (*out)[base + p];
        RawRow raw;
        raw.total = nl.total;
        raw.kept.reserve(nl.labels.size());
        for (LabelId id : nl.labels) {
            raw.kept.emplace_back(dict_[id].column,
                                  kind_ == LabelKind::HEADER ? dict_[id].seq_id : 0);
        }
        std::sort(raw.kept.begin(), raw.kept.end());
        cache_raw(keys[p], std::move(raw), (*costs)[base + p]);
    }
    return true;
}

void LabelRecorder::warm(const std::vector<node_index> &keys, DecodeBudget &budget,
                         ReadPacing *pacing) {
    if (costs_.size() != cache_.size())
        clear_cache();
    // as LabelQuery's: a cache of capacity zero keeps nothing, and its runs would not advance
    if (!max_cache_size_)
        return;
    // the misses, charged like the rest of the lookahead's read
    const uint64_t at_entry = budget.held();
    if (!budget.charge(buffer_bytes(keys.size(), sizeof(node_index))))
        return;
    std::vector<node_index> missing;
    missing.reserve(keys.size());
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key))
            missing.push_back(key);
    }
    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    DecodePacer &pacer = oracle_.pacer();
    const bool paced = pacing && pacer.target_ms > 0;
    size_t previous = 0;
    double previous_ms = 0;
    const size_t run = std::min(kMaxDecodeRun, max_cache_size_);
    for (size_t begin = 0; begin < missing.size(); begin += run) {
        const size_t len = std::min(run, missing.size() - begin);
        const uint64_t at_run = budget.held();
        if (!budget.charge(run_bytes_of(len, sizeof(RawRow))))
            break;
        std::vector<RawRow> raws(len);
        std::vector<uint64_t> raws_held(len);
        std::vector<KeyCost> run_costs(len);
        // as LabelQuery::warm: the run decoded in paced pieces and cached as one run
        bool ok = true;
        for (size_t at = 0; ok && at < len; ) {
            size_t piece = len - at;
            if (paced) {
                if (pacing->stop && pacing->stop()) {
                    pacing->interrupted = true;
                    budget.restore(at_entry);
                    return;
                }
                piece = pacer.next(piece, pacing->ms_left ? pacing->ms_left()
                                                          : std::numeric_limits<double>::infinity(),
                                   previous, previous_ms);
            }
            size_t built = 0;
            Timer timer;
            ok = decode_run(missing.data() + begin + at, piece, budget, &raws, &raws_held,
                            &run_costs, &built, at) == DecodeStatus::OK;
            previous_ms = timer.elapsed() * 1000;
            pacer.record(piece, previous_ms);
            previous = piece;
            at += piece;
        }
        if (!ok)
            break;      // the lookahead gives up within what is left: nothing depends on it
        uint64_t bytes = 0, kept = 0;
        for (const RawRow &raw : raws) {
            bytes += kCacheEntryBytes + kCostEntryBytes + raw_bytes(raw);
            kept += std::max<size_t>(1, raw.kept.size());
        }
        if (cache_bytes_ + bytes > max_cache_bytes_ || cache_.size() + len > max_cache_size_
                || cached_keys_ + kept > max_cache_keys_) {
            clear_cache();
        }
        if (bytes <= max_cache_bytes_ && len <= max_cache_size_ && kept <= max_cache_keys_) {
            for (size_t i = 0; i < len; ++i) {
                cache_raw(missing[begin + i], std::move(raws[i]), run_costs[i]);
            }
        }
        budget.restore(at_run);
    }
    budget.restore(at_entry);
}

} // namespace traversal
} // namespace graph
} // namespace mtg
