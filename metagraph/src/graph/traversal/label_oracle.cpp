#include "label_oracle.hpp"

#include <algorithm>
#include <set>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"


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

    // A CanonicalDBG wrapper means a PRIMARY index. Otherwise the mode tells native
    // CANONICAL graphs from BASIC ones.
    if (const auto *wrapper = dynamic_cast<const CanonicalDBG*>(&loaded)) {
        regime_ = Regime::PRIMARY;
        // a wrapper of our own: same node ids, own caches; the shared one is never mutated
        local_canonical_ = std::make_shared<CanonicalDBG>(wrapper->get_graph_ptr());
        if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&wrapper->get_graph()))
            local_canonical_->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));
        graph_ = local_canonical_.get();
    } else {
        if (loaded.get_mode() == DeBruijnGraph::PRIMARY) {
            throw std::invalid_argument("PRIMARY graphs must be wrapped into CanonicalDBG "
                                        "before traversal");
        }
        regime_ = loaded.get_mode() == DeBruijnGraph::CANONICAL ? Regime::CANONICAL
                                                               : Regime::BASIC;
        graph_ = &loaded;
    }
    key_mapper_ = std::make_unique<AnnotationKeyMapper>(*graph_);

    // get_matrix() may lock and flush (ColumnCompressed): call it once per request
    matrix_ = &anno_graph_.get_annotator().get_matrix();
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
    return node == npos ? npos : key_mapper_->key(node, kmer);
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

std::pair<uint64_t, uint64_t> LabelOracle::map_coord(Column column, uint64_t coord) const {
    assert(coord_to_header_);
    counters_.coords_mapped++;
    return coord_to_header_->map_single_coord(column, coord);
}

std::pair<uint64_t, uint64_t> LabelOracle::sequence_coords(Column column, uint64_t seq_id) const {
    assert(coord_to_header_);
    const uint64_t last = coord_to_header_->last_coord(column, seq_id);
    const uint64_t first = seq_id ? coord_to_header_->last_coord(column, seq_id - 1) + 1 : 0;
    return { first, last };
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


LabelQuery::LabelQuery(const LabelOracle &oracle,
                       std::vector<LabelRef> labels,
                       size_t max_cache_size)
      : oracle_(oracle), labels_(std::move(labels)), max_cache_size_(max_cache_size) {
    if (labels_.empty())
        throw std::invalid_argument("Empty permitted label set");

    std::set<std::pair<Column, uint64_t>> headers;
    for (LabelId i = 0; i < labels_.size(); ++i) {
        const LabelRef &ref = labels_[i];
        if (ref.kind == LabelKind::COLUMN) {
            if (!column_labels_.try_emplace(ref.column, i).second)
                throw std::invalid_argument("Duplicate label '" + ref.name + "'");
            continue;
        }
        if (!oracle_.coord_to_header())
            throw std::invalid_argument("Sequence header labels require an annotation with "
                                        "coordinates and a CoordToHeader index");
        if (!headers.emplace(ref.column, ref.seq_id).second)
            throw std::invalid_argument("Duplicate label '" + ref.name + "'");
        const auto [first, last] = oracle_.sequence_coords(ref.column, ref.seq_id);
        header_ranges_[ref.column].push_back(HeaderRange{ first, last, i });
    }
    for (auto it = header_ranges_.begin(); it != header_ranges_.end(); ++it) {
        std::sort(it.value().begin(), it.value().end(),
                  [](const HeaderRange &a, const HeaderRange &b) { return a.first < b.first; });
    }
    tuples_ = !header_ranges_.empty();
}

void LabelQuery::hits_from_row(const BinaryMatrix::SetBitPositions &row, NodeHits *hits) const {
    for (Column c : row) {
        auto it = column_labels_.find(c);
        if (it != column_labels_.end())
            hits->push_back(it->second);
    }
    std::sort(hits->begin(), hits->end());
}

void LabelQuery::hits_from_tuples(const MultiIntMatrix::RowTuples &row, NodeHits *hits) const {
    for (const auto &[c, coords] : row) {
        if (auto it = column_labels_.find(c); it != column_labels_.end())
            hits->push_back(it->second);
        auto rt = header_ranges_.find(c);
        if (rt == header_ranges_.end())
            continue;
        const std::vector<HeaderRange> &ranges = rt->second;
        for (uint64_t coord : coords) {
            // the last range starting at or before the coordinate (the ranges are disjoint)
            auto it = std::upper_bound(ranges.begin(), ranges.end(), coord,
                                       [](uint64_t x, const HeaderRange &r) { return x < r.first; });
            if (it != ranges.begin() && coord <= (--it)->last)
                hits->push_back(it->label);
        }
    }
    std::sort(hits->begin(), hits->end());
    hits->erase(std::unique(hits->begin(), hits->end()), hits->end());
}

void LabelQuery::fetch_uncached(const std::vector<node_index> &keys) {
    if (keys.empty())
        return;
    if (cache_.size() + keys.size() > max_cache_size_)
        cache_.clear();

    std::vector<Row> rows;
    rows.reserve(keys.size());
    for (node_index key : keys) {
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    if (tuples_) {
        auto tuple_rows = oracle_.get_row_tuples(rows);
        for (size_t i = 0; i < keys.size(); ++i) {
            hits_from_tuples(tuple_rows[i], &cache_[keys[i]]);
        }
    } else {
        auto set_rows = oracle_.get_rows(rows);
        for (size_t i = 0; i < keys.size(); ++i) {
            hits_from_row(set_rows[i], &cache_[keys[i]]);
        }
    }
}

std::vector<LabelQuery::NodeHits> LabelQuery::fetch(const std::vector<node_index> &keys) {
    auto &counters = oracle_.counters();
    std::vector<node_index> missing;
    tsl::hopscotch_set<node_index> requested;
    for (node_index key : keys) {
        if (key == npos)
            continue;
        counters.rows_requested++;
        if (cache_.count(key)) {
            counters.cache_hits++;
        } else if (requested.insert(key).second) {
            missing.push_back(key);
        }
    }
    // a clear in fetch_uncached may drop keys requested here: read them all again then
    if (cache_.size() + missing.size() > max_cache_size_) {
        missing.clear();
        requested.clear();
        for (node_index key : keys) {
            if (key != npos && requested.insert(key).second)
                missing.push_back(key);
        }
    }
    fetch_uncached(missing);

    std::vector<NodeHits> result(keys.size());
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] != npos)
            result[i] = cache_.at(keys[i]);
    }
    return result;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
