#include "label_oracle.hpp"

#include <algorithm>
#include <mutex>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;
using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;


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

namespace {

// header -> (column, seq_id), built once per loaded CoordToHeader
struct HeaderIndex {
    tsl::hopscotch_map<std::string, std::pair<Column, uint64_t>> map;
};

std::shared_ptr<const HeaderIndex> get_header_index(const annot::CoordToHeader &cth) {
    static std::mutex mu;
    static std::vector<std::pair<const annot::CoordToHeader*, std::shared_ptr<const HeaderIndex>>> cache;
    std::lock_guard<std::mutex> lock(mu);
    for (const auto &[ptr, index] : cache) {
        if (ptr == &cth)
            return index;
    }
    auto index = std::make_shared<HeaderIndex>();
    size_t total = 0;
    for (Column col = 0; col < cth.num_columns(); ++col) {
        total += cth.num_sequences(col);
    }
    index->map.reserve(total);
    for (Column col = 0; col < cth.num_columns(); ++col) {
        const auto &headers = cth.get_headers(col);
        for (uint64_t s = 0; s < headers.size(); ++s) {
            index->map.try_emplace(headers[s], col, s);
        }
    }
    logger->trace("Built header index with {} sequences in {} columns", total, cth.num_columns());
    cache.emplace_back(&cth, index);
    return index;
}

} // namespace

std::optional<LabelRef> LabelOracle::find_header(const std::string &name) const {
    if (!coord_to_header_)
        return std::nullopt;
    auto index = get_header_index(*coord_to_header_);
    auto it = index->map.find(name);
    if (it == index->map.end())
        return std::nullopt;
    LabelRef ref;
    ref.kind = LabelKind::HEADER;
    ref.column = it->second.first;
    ref.seq_id = it->second.second;
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

std::vector<BinaryMatrix::SetBitPositions>
LabelOracle::get_rows(const std::vector<Row> &rows) const {
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
    Timer timer;
    auto result = tuples_->get_row_tuples(rows);
    counters_.tuple_rows_fetched += rows.size();
    counters_.fetch_seconds += timer.elapsed();
    return result;
}

bool LabelOracle::get(Row row, Column column) const {
    assert(get_entry_);
    counters_.direct_reads++;
    return get_entry_->get(row, column);
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
            if (!header_labels_[ref.column].try_emplace(ref.seq_id, i).second)
                throw std::invalid_argument("Duplicate label '" + ref.name + "'");
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
    for (const auto &[c, coords] : row) {
        if (auto it = column_labels_.find(c); it != column_labels_.end()) {
            Hit hit{ it->second, {} };
            if (with_coords_) {
                hit.coords.assign(coords.begin(), coords.end());
                std::sort(hit.coords.begin(), hit.coords.end());
            }
            hits->push_back(std::move(hit));
        }
        auto ht = header_labels_.find(c);
        if (ht == header_labels_.end())
            continue;
        for (Coord coord : coords) {
            auto [seq_id, local] = oracle_.map_coord(c, coord);
            auto lt = ht->second.find(seq_id);
            if (lt == ht->second.end())
                continue;
            auto found = std::find_if(hits->begin(), hits->end(),
                                      [&](const Hit &h) { return h.label == lt->second; });
            if (found == hits->end()) {
                hits->push_back(Hit{ lt->second, {} });
                found = hits->end() - 1;
            }
            if (with_coords_)
                found->coords.push_back(local);
        }
    }
    for (Hit &hit : *hits) {
        std::sort(hit.coords.begin(), hit.coords.end());
    }
    std::sort(hits->begin(), hits->end(),
              [](const Hit &a, const Hit &b) { return a.label < b.label; });
}

void LabelQuery::fetch_uncached(const std::vector<node_index> &keys) {
    if (keys.empty())
        return;

    std::vector<Row> rows;
    rows.reserve(keys.size());
    for (node_index key : keys) {
        assert(key != npos);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }

    std::vector<NodeHits> result(keys.size());
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
            break;
        }
        case Path::TUPLES: {
            auto fetched = oracle_.get_row_tuples(rows);
            for (size_t i = 0; i < rows.size(); ++i) {
                hits_from_tuples(fetched[i], &result[i]);
            }
            break;
        }
    }

    // Always insert: a silent no-op here would make the caller's cache_.at() throw.
    // The caller evicts before fetching when the cache would overflow, so the cache can
    // exceed max_cache_size_ only by the working set of a single call.
    for (size_t i = 0; i < keys.size(); ++i) {
        cache_[keys[i]] = std::move(result[i]);
    }
}

std::vector<LabelQuery::NodeHits> LabelQuery::fetch(const std::vector<node_index> &keys) {
    // Distinct keys this call must be able to answer. Eviction below may drop entries
    // that were cached on entry, so the whole working set — not just the misses — has
    // to be (re)fetched; otherwise the lookups at the end would throw.
    std::vector<node_index> wanted;
    wanted.reserve(keys.size());
    for (node_index key : keys) {
        if (key == npos)
            continue;
        oracle_.counters().rows_requested++;
        if (cache_.count(key))
            oracle_.counters().cache_hits++;
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
        if (cache_.size() + missing.size() > max_cache_size_) {
            // a walk moves forward, so evict wholesale rather than tracking recency,
            // then refetch everything this call needs
            cache_.clear();
            missing = wanted;
        }
        fetch_uncached(missing);
    }

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

void LabelQuery::warm(const std::vector<node_index> &keys) {
    std::vector<node_index> missing;
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key))
            missing.push_back(key);
    }
    if (missing.empty())
        return;
    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    if (cache_.size() + missing.size() > max_cache_size_)
        cache_.clear();
    if (missing.size() >= max_cache_size_)
        return;     // would not fit even into an empty cache: nothing to warm
    fetch_uncached(missing);
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
    if (cache_.size() >= max_cache_size_)
        cache_.clear();
    fetch_uncached({ key });
    return cache_.at(key);
}

} // namespace traversal
} // namespace graph
} // namespace mtg
