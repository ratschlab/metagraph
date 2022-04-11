#include "annotation_buffer.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"

namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::binmat::BinaryMatrix::Row Row;
typedef annot::binmat::BinaryMatrix::Column Column;

// dummy index for an unfetched annotations
static constexpr size_t nannot = std::numeric_limits<size_t>::max();

template <class AIt, class BIt, class CIt, class OutIt, class OutIt2>
void matched_intersection(AIt a_begin, AIt a_end, BIt a_c_begin,
                          CIt b_begin, CIt b_end,
                          OutIt out1, OutIt2 out2) {
    while (a_begin != a_end && b_begin != b_end) {
        if (*a_begin < *b_begin) {
            ++a_begin;
            ++a_c_begin;
        } else if (*a_begin > *b_begin) {
            ++b_begin;
        } else {
            *out1 = *a_begin;
            *out2 = *a_c_begin;
            ++out1;
            ++out2;
            ++a_begin;
            ++a_c_begin;
            ++b_begin;
        }
    }
}

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph,
                                   const Annotator &annotator,
                                   size_t row_batch_size,
                                   size_t max_coords_buffered)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        canonical_(dynamic_cast<const CanonicalDBG*>(&graph_)),
        column_sets_({ {} }),
        label_coords_cache_(max_coords_buffered),
        row_batch_size_(row_batch_size) {
    if (multi_int_ && graph_.get_mode() != DeBruijnGraph::BASIC) {
        multi_int_ = nullptr;
        logger->warn("Coordinates not supported when aligning to CANONICAL "
                     "or PRIMARY mode graphs");
    }
}

#ifndef NDEBUG
bool check_coords(const annot::matrix::MultiIntMatrix &mat,
                  annot::binmat::BinaryMatrix::Row row,
                  const AnnotationBuffer::Columns &columns,
                  const std::optional<const AnnotationBuffer::CoordinateSet> &coord_set = std::nullopt) {
    assert(!coord_set || columns.size() == coord_set->size());
    auto row_tuples = mat.get_row_tuples(row);
    assert(row_tuples.size() == columns.size());

    std::sort(row_tuples.begin(), row_tuples.end());
    for (size_t i = 0; i < row_tuples.size(); ++i) {
        const auto &[c, tuple] = row_tuples[i];
        assert(c == columns[i]);
        assert(std::is_sorted(tuple.begin(), tuple.end()));
        assert(!coord_set
            || std::equal(tuple.begin(), tuple.end(), (*coord_set)[i].begin(), (*coord_set)[i].end()));
    }

    return true;
}

template <class T>
inline T sorted(T&& v) {
    std::sort(v.begin(), v.end());
    return v;
}

#endif

auto AnnotationBuffer
::get_label_and_coord_diffs(node_index node, const std::vector<node_index> &nexts, bool flipped)
        -> std::vector<std::pair<Columns, std::shared_ptr<CoordinateSet>>> {
    node_index node_base = node;
    if (canonical_) {
        node_base = canonical_->get_base_node(node);
    } else if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        node_base = map_to_nodes(graph_, graph_.get_node_sequence(node))[0];
    }

    std::vector<std::pair<Columns, std::shared_ptr<CoordinateSet>>> ret_vals;
    ret_vals.reserve(nexts.size());

    for (node_index next : nexts) {
        node_index a = node_base;
        auto find_a = node_to_cols_.find(a);
        assert(find_a != node_to_cols_.end());
        size_t labels_a = find_a->second;

        node_index b = next;
        if (canonical_) {
            b = canonical_->get_base_node(next);
        } else if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
            b = map_to_nodes(graph_, graph_.get_node_sequence(next))[0];
        }

        if (flipped || (a != node && b != next)) {
            std::swap(a, b);
            flipped = true;
        }

        auto row_a = AnnotatedDBG::graph_to_anno_index(a);
        auto row_b = AnnotatedDBG::graph_to_anno_index(b);

        auto &[columns, coords] = ret_vals.emplace_back();

        if (has_coordinates()) {
            Columns next_columns;
            CoordinateSet next_coords;
            auto diff = multi_int_->get_row_tuple_diff(row_a, row_b);
            assert(std::is_sorted(diff.begin(), diff.end()));
            next_columns.reserve(diff.size());
            next_coords.reserve(diff.size());
            for (auto&& [c, tuple] : diff) {
                assert(std::is_sorted(tuple.begin(), tuple.end()));
                next_columns.push_back(c);
                next_coords.emplace_back(tuple.begin(), tuple.end());
            }

            std::swap(next_columns, columns);
            coords = std::make_shared<CoordinateSet>(std::move(next_coords));
        } else {
            columns = annotator_.get_matrix().get_diff(row_a, row_b);
            assert(std::is_sorted(columns.begin(), columns.end()));
        }

        if (flipped) {
            std::swap(a, b);
#ifndef NDEBUG
            std::swap(row_a, row_b);
#endif
        }

        auto find_b = node_to_cols_.find(b);
        if (find_b == node_to_cols_.end() || find_b->second == nannot) {
            if (columns.size()) {
                const auto &node_columns = column_sets_.data()[labels_a];
                Columns next_columns;
                if (coords) {
                    auto node_coords = get_coords_from_it(find_a, false);
                    if (!flipped) {
                        for (auto &tuple : *node_coords) {
                            for (auto &c : tuple) {
                                ++c;
                            }
                        }
                    }
                    utils::match_indexed_values(
                        node_columns.begin(), node_columns.end(), node_coords->begin(),
                        columns.begin(), columns.end(), coords->begin(),
                        [&](Column c, const auto &coords, const auto &other_coords) {
                            // if other_coords is empty, then it means that there
                            // was a column in next, but not in node
                            assert(coords.size());
                            assert(flipped || other_coords.size());
                            if (other_coords.size() && coords != other_coords)
                                next_columns.push_back(c);
                        },
                        [&](Column c, const auto &) { next_columns.push_back(c); },
                        [&](Column c, const auto &) { next_columns.push_back(c); }
                    );
                    if (flipped) {
                        for (auto &tuple : *coords) {
                            for (auto &c : tuple) {
                                ++c;
                            }
                        }
                    }
                } else {
                    std::set_symmetric_difference(node_columns.begin(), node_columns.end(),
                                                  columns.begin(), columns.end(),
                                                  std::back_inserter(next_columns));
                }

                assert(next_columns == sorted(annotator_.get_matrix().get_row(row_b)));

                labels_a = next_columns.size()
                    ? cache_column_set(std::move(next_columns))
                    : 0;
            }


            auto [it, inserted] = node_to_cols_.try_emplace(b, labels_a);
            if (!inserted)
                it.value() = labels_a;

            if (b != next) {
                auto [it, inserted] = node_to_cols_.try_emplace(next, labels_a);
                if (!inserted)
                    it.value() = labels_a;
            }
        }
    }
    return ret_vals;
}

void AnnotationBuffer::fetch_queued_annotations() {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    const DeBruijnGraph *base_graph = &graph_;

    if (canonical_)
        base_graph = &canonical_->get_graph();

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    auto fetch_row_batch = [&](auto&& queued_nodes, auto&& queued_rows) {
        if (queued_nodes.empty())
            return;

        auto node_it = queued_nodes.begin();
        auto row_it = queued_rows.begin();
        if (has_coordinates()) {
            for (auto&& tuples : multi_int_->get_row_tuples(queued_rows)) {
                std::sort(tuples.begin(), tuples.end());
                std::all_of(tuples.begin(), tuples.end(), [&](const auto &a) {
                    return a.second.size();
                });
                assert(node_it != queued_nodes.end());
                assert(node_to_cols_.count(*node_it));
                assert(node_to_cols_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));

                Columns columns;
                CoordinateSet coord_set;
                columns.reserve(tuples.size());
                coord_set.reserve(tuples.size());
                for (auto&& [c, tuple] : tuples) {
                    columns.push_back(c);
                    coord_set.emplace_back(tuple.begin(), tuple.end());
                }
                size_t label_i = cache_column_set(std::move(columns));
                assert(AnnotatedDBG::anno_to_graph_index(*row_it) == *node_it
                        && "coordinates only supported for BASIC graphs");
                auto find = node_to_cols_.find(*node_it);
                assert(find != node_to_cols_.end());
                find.value() = label_i;
                label_coords_cache_.Put(find - node_to_cols_.begin(), std::move(coord_set));
                ++node_it;
                ++row_it;
            }
        } else {
            for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
                std::sort(labels.begin(), labels.end());
                assert(node_it != queued_nodes.end());
                assert(node_to_cols_.count(*node_it));
                assert(node_to_cols_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));

                size_t label_i = cache_column_set(std::move(labels));
                node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
                if (canonical_) {
                    node_to_cols_[base_node] = label_i;
                } else {
                    node_to_cols_[*node_it] = label_i;
                    if (base_node != *node_it) {
                        assert(graph_.get_mode() != DeBruijnGraph::BASIC);
                        node_to_cols_[base_node] = label_i;
                    }
                }
                ++node_it;
                ++row_it;
            }
        }
    };

    std::vector<node_index> queued_nodes;
    std::vector<Row> queued_rows;
    for (const auto &path : queued_paths_) {
        std::vector<node_index> base_path;
        if (base_graph->get_mode() == DeBruijnGraph::CANONICAL) {
            // TODO: avoid this call of spell_path
            std::string query = spell_path(graph_, path);
            base_path = map_to_nodes(*base_graph, query);

        } else if (canonical_) {
            base_path.reserve(path.size());
            for (node_index node : path) {
                base_path.emplace_back(canonical_->get_base_node(node));
            }

        } else {
            assert(graph_.get_mode() == DeBruijnGraph::BASIC);
            base_path = path;
            if (dynamic_cast<const RCDBG*>(&graph_))
                std::reverse(base_path.begin(), base_path.end());
        }

        assert(base_path.size() == path.size());

        for (size_t i = 0; i < path.size(); ++i) {
            if (base_path[i] == DeBruijnGraph::npos) {
                // this can happen when the base graph is CANONICAL and path[i] is a
                // dummy node
                node_to_cols_.try_emplace(path[i], 0);
                continue;
            }

            if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i]))) {
                // skip dummy nodes
                node_to_cols_.try_emplace(base_path[i], 0);
                if (graph_.get_mode() == DeBruijnGraph::CANONICAL && base_path[i] != path[i])
                    node_to_cols_.try_emplace(path[i], 0);

                continue;
            }

            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            if (canonical_ || graph_.get_mode() == DeBruijnGraph::BASIC) {
                if (node_to_cols_.try_emplace(base_path[i], nannot).second) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                    if (queued_rows.size() >= row_batch_size_) {
                        fetch_row_batch(std::move(queued_nodes), std::move(queued_rows));
                        queued_nodes = decltype(queued_nodes){};
                        queued_rows = decltype(queued_rows){};
                    }
                }

                continue;
            }

            assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);

            auto find_a = node_to_cols_.find(path[i]);
            auto find_b = node_to_cols_.find(base_path[i]);

            if (find_a == node_to_cols_.end() && find_b == node_to_cols_.end()) {
                node_to_cols_.try_emplace(path[i], nannot);
                queued_rows.push_back(row);
                queued_nodes.push_back(path[i]);

                if (path[i] != base_path[i]) {
                    node_to_cols_.emplace(base_path[i], nannot);
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                }
            } else if (find_a == node_to_cols_.end() && find_b != node_to_cols_.end()) {
                node_to_cols_.try_emplace(path[i], find_b->second);
                if (find_b->second == nannot) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(path[i]);
                }
            } else if (find_a != node_to_cols_.end() && find_b == node_to_cols_.end()) {
                node_to_cols_.try_emplace(base_path[i], find_a->second);
            } else {
                size_t label_i = std::min(find_a->second, find_b->second);
                if (label_i != nannot) {
                    find_a.value() = label_i;
                    find_b.value() = label_i;
                }
            }

            if (queued_rows.size() >= row_batch_size_) {
                fetch_row_batch(std::move(queued_nodes), std::move(queued_rows));
                queued_nodes = decltype(queued_nodes){};
                queued_rows = decltype(queued_rows){};
            }
        }
    }

    fetch_row_batch(std::move(queued_nodes), std::move(queued_rows));

    queued_paths_.clear();
}

void AnnotationBuffer::prefetch_coords(const std::vector<node_index> &nodes) const {
    if (!has_coordinates() || nodes.empty())
        return;

    assert(graph_.get_mode() == DeBruijnGraph::BASIC && "coordinates only supported for BASIC graphs");

    std::vector<Row> rows;
    rows.reserve(nodes.size());
    std::vector<size_t> indices;
    indices.reserve(nodes.size());
    std::vector<size_t> label_is;
    label_is.reserve(nodes.size());
    for (node_index node : nodes) {
        rows.emplace_back(AnnotatedDBG::graph_to_anno_index(node));

        auto it = node_to_cols_.find(node);
        assert(it != node_to_cols_.end());
        assert(it->second != nannot);
        indices.emplace_back(it - node_to_cols_.begin());
        label_is.emplace_back(it->second);
    }

    if (std::all_of(indices.begin(), indices.end(), [&](size_t index) {
                return label_coords_cache_.Cached(index);
            })) {
        return;
    }

    auto id_it = indices.begin();
    auto label_it = label_is.begin();
    const Columns *columns = &column_sets_.data()[*label_it];

    auto fetch = label_coords_cache_.TryGet(*id_it);
    assert(check_coords(*multi_int_, rows[0], *columns, fetch));
    if (!fetch) {
        auto first_tuple = multi_int_->get_row_tuples(rows[0]);
        std::all_of(first_tuple.begin(), first_tuple.end(), [&](const auto &a) {
            return a.second.size();
        });
        std::sort(first_tuple.begin(), first_tuple.end());
        CoordinateSet coord_set;
        coord_set.reserve(first_tuple.size());
        for (auto&& [column, coords] : first_tuple) {
            coord_set.emplace_back(coords.begin(), coords.end());
        }
        fetch = coord_set;
        assert(check_coords(*multi_int_, rows[0], *columns, fetch));
        label_coords_cache_.Put(*id_it, std::move(coord_set));
    }

    // fetch diffs instead
    assert(fetch);
    assert(columns->size() == fetch->size());
    annot::matrix::MultiIntMatrix::RowTuples first_tuples;
    first_tuples.reserve(columns->size());
    for (size_t i = 0; i < columns->size(); ++i) {
        first_tuples.emplace_back((*columns)[i], annot::matrix::MultiIntMatrix::Tuple((*fetch)[i].begin(), (*fetch)[i].end()));
    }

    ++id_it;
    ++label_it;
    auto tuple_diffs = multi_int_->get_row_tuple_diffs(rows, &first_tuples);
    assert(tuple_diffs[0] == first_tuples);

    for (size_t i = 1; i < tuple_diffs.size(); ++i) {
        auto &tuple_diff = tuple_diffs[i];
        assert(std::is_sorted(tuple_diff.begin(), tuple_diff.end()));
        columns = &column_sets_.data()[*label_it];
        assert(*columns == multi_int_->get_row(nodes[i]));

        for (auto &tuple : *fetch) {
            for (auto &c : tuple) {
                ++c;
            }
        }
        if (tuple_diff.size()) {
            CoordinateSet next_coord_set;
            auto it = columns->begin();
            auto it_c = fetch->begin();
            auto jt = tuple_diff.begin();
            while (it != columns->end() && jt != tuple_diff.end()) {
                if (*it < jt->first) {
                    next_coord_set.push_back(*it_c);
                    ++it;
                    ++it_c;
                } else if (*it > jt->first) {
                    next_coord_set.emplace_back(jt->second.begin(), jt->second.end());
                    ++jt;
                } else {
                    next_coord_set.emplace_back();
                    auto c_it = it_c->begin();
                    auto d_it = jt->second.begin();
                    while (c_it != it_c->end() && d_it != jt->second.end()) {
                        int64_t d_coord = *d_it;
                        if (*c_it < d_coord) {
                            next_coord_set.back().emplace_back(*c_it);
                            ++c_it;
                        } else if (*c_it > d_coord) {
                            next_coord_set.back().emplace_back(d_coord);
                            ++d_it;
                        }
                    }
                    std::copy(c_it, it_c->end(), std::back_inserter(next_coord_set.back()));
                    std::copy(d_it, jt->second.end(), std::back_inserter(next_coord_set.back()));
                    ++it;
                    ++it_c;
                    ++jt;
                }

            }
            std::copy(it_c, fetch->end(), std::back_inserter(next_coord_set));
            std::transform(jt, tuple_diff.end(), std::back_inserter(next_coord_set),
                           [&](const auto &t) { return CoordinateSet::value_type(t.second.begin(), t.second.end()); });
            std::swap(*fetch, next_coord_set);
        }
        assert(check_coords(*multi_int_, rows[i], *columns, *fetch));
        label_coords_cache_.Put(*id_it, *fetch);
        ++id_it;
        ++label_it;
    }
}

auto AnnotationBuffer::get_labels_and_coords(node_index node, bool skip_unfetched) const
        -> std::pair<const Columns*, std::shared_ptr<CoordinateSet>> {
    if (canonical_)
        node = canonical_->get_base_node(node);

    auto it = node_to_cols_.find(node);

    // if the node hasn't been seen before, or if its annotations haven't
    // been fetched, return nothing
    if (it == node_to_cols_.end() || it->second == nannot)
        return {};

    return std::make_pair(&column_sets_.data()[it->second],
                          has_coordinates() ? get_coords_from_it(it, skip_unfetched)
                                            : std::shared_ptr<CoordinateSet>{});
}

auto AnnotationBuffer::get_coords_from_it(VectorMap<node_index, size_t>::const_iterator it,
                                          bool skip_unfetched,
                                          const Columns *label_subset) const
        -> std::shared_ptr<CoordinateSet> {
    assert(has_coordinates());
    assert(it != node_to_cols_.end());
    assert(it->second != nannot);

    const auto &columns = column_sets_.data()[it->second];

    if (columns.empty())
        return {};

    size_t index = it - node_to_cols_.begin();

    auto coord_set = label_coords_cache_.TryGet(index);
    if (!coord_set) {
        if (skip_unfetched)
            return {};

        auto fetch
            = multi_int_->get_row_tuples(AnnotatedDBG::graph_to_anno_index(it->first));
        std::sort(fetch.begin(), fetch.end());
        coord_set = CoordinateSet{};
        assert(fetch.size() == columns.size());
        for (size_t i = 0; i < fetch.size(); ++i) {
            auto &[column, coords] = fetch[i];
            assert(column == columns[i]);
            coord_set->emplace_back(coords.begin(), coords.end());
        }

        label_coords_cache_.Put(index, *coord_set);
    }

    assert(check_coords(*multi_int_, AnnotatedDBG::graph_to_anno_index(it->first), columns, *coord_set));

    if (!label_subset)
        return std::make_shared<CoordinateSet>(std::move(*coord_set));

    Columns dummy;
    CoordinateSet coord_subset;
    matched_intersection(columns.begin(), columns.end(), coord_set->begin(),
                         label_subset->begin(), label_subset->end(),
                         std::back_inserter(dummy),
                         std::back_inserter(coord_subset));
    return std::make_shared<CoordinateSet>(std::move(coord_subset));
}

auto AnnotationBuffer::get_coords(node_index node,
                                  bool skip_unfetched,
                                  const Columns *label_subset) const
        -> std::shared_ptr<CoordinateSet> {
    if (!has_coordinates())
        return {};

    assert(graph_.get_mode() == DeBruijnGraph::BASIC && "coordinates only supported for BASIC graphs");

    auto it = node_to_cols_.find(node);

    // if the node hasn't been seen before, or if its annotations haven't
    // been fetched, return nothing
    if (it == node_to_cols_.end() || it->second == nannot)
        return {};

    return get_coords_from_it(it, skip_unfetched, label_subset);
}

} // namespace align
} // namespace graph
} // namespace mtg
