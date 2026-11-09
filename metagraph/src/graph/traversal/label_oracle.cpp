#include "label_oracle.hpp"

#include <algorithm>
#include <set>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;
using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;


LabelOracle::LabelOracle(const AnnotatedDBG &anno_graph) : anno_graph_(anno_graph) {
    const DeBruijnGraph &loaded = anno_graph_.get_graph();

    // A CanonicalDBG wrapper means a PRIMARY index. Otherwise the mode tells native
    // CANONICAL graphs from BASIC ones.
    if (const auto *wrapper = dynamic_cast<const CanonicalDBG*>(&loaded)) {
        regime_ = Regime::PRIMARY;
        // a wrapper of our own: same node ids, own caches; the shared one is never mutated
        local_canonical_ = CanonicalDBG::with_node_first_cache(wrapper->get_graph_ptr());
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

    // get_matrix() may lock and flush (ColumnCompressed): call it once per request
    matrix_ = &anno_graph_.get_annotator().get_matrix();
    tuples_ = dynamic_cast<const MultiIntMatrix*>(matrix_);
    row_diff_ = dynamic_cast<const annot::matrix::IRowDiff*>(matrix_);
    coord_to_header_ = anno_graph_.get_coord_to_header();
    if (coord_to_header_ && !tuples_) {
        logger->warn("CoordToHeader present but the annotation has no coordinates; "
                     "header labels are unavailable");
        coord_to_header_ = nullptr;
    }
}

LabelOracle::~LabelOracle() {}

size_t LabelOracle::get_k() const { return graph_->get_k(); }

uint64_t LabelOracle::num_columns() const { return matrix_->num_columns(); }

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

bool LabelOracle::has_header_records(Column column) const {
    return coord_to_header_ && column < coord_to_header_->num_columns()
            && coord_to_header_->num_sequences(column);
}

std::pair<uint64_t, uint64_t> LabelOracle::map_coord(Column column, uint64_t coord) const {
    assert(coord_to_header_);
    return coord_to_header_->map_single_coord(column, coord);
}


bool LabelOracle::is_ambiguous(const std::string &name) const {
    return find_header(name)
        && (find_column(name) || coord_to_header_->is_repeated_header(name));
}

LabelRef LabelOracle::resolve_label(const std::string &name) const {
    auto column = find_column(name);
    if (is_ambiguous(name)) {
        throw std::invalid_argument("Label '" + name + "' is ambiguous: " + (column
                ? "it names both an annotation column and a sequence header"
                : "several sequence records have this header"));
    }
    auto header = find_header(name);
    if (column) {
        LabelRef ref;
        ref.kind = LabelKind::COLUMN;
        ref.column = *column;
        ref.name = name;
        return ref;
    }
    if (header)
        return *header;
    throw std::invalid_argument("Label not found in the annotation: '" + name + "'"
                                + (coord_to_header_ ? "" : " (no sequence header index loaded)"));
}

std::vector<BinaryMatrix::SetBitPositions>
LabelOracle::get_rows(const std::vector<Row> &rows,
                      annot::matrix::RowDiffCache<BinaryMatrix::SetBitPositions> *cache) const {
    std::vector<BinaryMatrix::SetBitPositions> result;
    if (!cache || !row_diff_ || !row_diff_->get_rows_cached(rows, *cache, &result))
        result = matrix_->get_rows(rows);
    return result;
}

std::vector<MultiIntMatrix::RowTuples>
LabelOracle::get_row_tuples(const std::vector<Row> &rows,
                            annot::matrix::RowDiffCache<MultiIntMatrix::RowTuples> *cache) const {
    if (!tuples_)
        throw std::logic_error("The annotation has no coordinates");
    std::vector<MultiIntMatrix::RowTuples> result;
    if (!cache || !row_diff_ || !row_diff_->get_row_tuples_cached(rows, *cache, &result))
        result = tuples_->get_row_tuples(rows);
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
        const auto [first, last] = oracle_.coord_to_header()->sequence_coords(ref.column, ref.seq_id);
        header_ranges_[ref.column].push_back(HeaderRange{ first, last, i });
    }
    for (auto it = header_ranges_.begin(); it != header_ranges_.end(); ++it) {
        std::sort(it.value().begin(), it.value().end(),
                  [](const HeaderRange &a, const HeaderRange &b) { return a.first < b.first; });
    }
    // A row-diff coordinate annotation decodes tuple rows for rows too: reading them lets
    // the decoded rows be kept (RowDiffCache) for column labels as well.
    tuples_ = !header_ranges_.empty() || (oracle_.has_coordinates() && oracle_.is_row_diff());
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

void LabelQuery::read_hits(const std::vector<node_index> &keys,
                           tsl::hopscotch_map<node_index, NodeHits> *out) {
    if (keys.empty())
        return;
    std::vector<Row> rows;
    rows.reserve(keys.size());
    for (node_index key : keys) {
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    if (tuples_) {
        auto tuple_rows = oracle_.get_row_tuples(rows, &tuple_cache_);
        for (size_t i = 0; i < keys.size(); ++i) {
            hits_from_tuples(tuple_rows[i], &(*out)[keys[i]]);
        }
    } else {
        auto set_rows = oracle_.get_rows(rows, &row_cache_);
        for (size_t i = 0; i < keys.size(); ++i) {
            hits_from_row(set_rows[i], &(*out)[keys[i]]);
        }
    }
}

std::vector<LabelQuery::NodeHits> LabelQuery::fetch(const std::vector<node_index> &keys) {
    std::vector<node_index> missing;
    tsl::hopscotch_set<node_index> requested;
    for (node_index key : keys) {
        if (key != npos && !cache_.count(key) && requested.insert(key).second)
            missing.push_back(key);
    }
    // when the new keys do not fit, clear the cache and read every key of this batch
    if (cache_.size() + missing.size() > max_cache_size_) {
        cache_.clear();
        missing.clear();
        requested.clear();
        for (node_index key : keys) {
            if (key != npos && requested.insert(key).second)
                missing.push_back(key);
        }
    }

    // a batch larger than the cache is answered without caching it
    tsl::hopscotch_map<node_index, NodeHits> uncached;
    auto *found = missing.size() > max_cache_size_ ? &uncached : &cache_;
    read_hits(missing, found);

    std::vector<NodeHits> result(keys.size());
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] != npos)
            result[i] = found->at(keys[i]);
    }
    return result;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
