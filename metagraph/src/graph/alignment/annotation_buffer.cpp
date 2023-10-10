#include "annotation_buffer.hpp"

#include <tsl/hopscotch_set.h>

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_set.hpp"
#include "common/algorithms.hpp"
#include "graph/graph_extensions/hll_wrapper.hpp"

namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::matrix::BinaryMatrix::Row Row;
typedef annot::matrix::BinaryMatrix::Column Column;

// dummy index for an unfetched annotations
static constexpr size_t nannot = std::numeric_limits<size_t>::max();

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        canonical_(dynamic_cast<const CanonicalDBG*>(&graph_)),
        column_sets_({ {} }) {
    if (multi_int_ && graph_.get_mode() != DeBruijnGraph::BASIC) {
        logger->warn("Coordinates not supported when aligning to CANONICAL "
                     "or PRIMARY mode graphs. Alignments may be truncated.");
    }
}

bool AnnotationBuffer::labels_valid(const Alignment &alignment) const {
    for (size_t i = 0; i < alignment.get_nodes().size(); ++i) {
        const auto &labels = alignment.get_columns(i);
        if (!check_node_labels_is_superset(labels, { alignment.get_nodes()[i] }))
            return false;
    }

    return true;
}

bool AnnotationBuffer
::check_node_labels_is_superset(const Columns &c, const std::vector<node_index> &nodes) const {
    if (c.empty())
        return true;

    for (node_index node : nodes) {
        const auto *labels = get_labels(node);
        if (!labels) {
            logger->error("Labels for node {} ({}) have not been fetched",
                          node, canonical_ ? canonical_->get_base_node(node) : node);
            return false;
        }

        Columns diff;
        std::set_difference(c.begin(), c.end(), labels->begin(), labels->end(),
                            std::back_inserter(diff));
        if (diff.size()) {
            std::vector<std::string> diff_labels;
            diff_labels.reserve(diff.size());
            const auto &label_encoder = annotator_.get_label_encoder();
            for (auto c : diff) {
                diff_labels.emplace_back(label_encoder.decode(c));
            }
            logger->error("Node {} does not have labels: {}", node, fmt::join(diff_labels, ";"));
            return false;
        }
    }

    return true;
}

void AnnotationBuffer::fetch_annotations(const std::vector<std::vector<node_index>> &queued_paths,
                                         VectorSet<Columns, utils::VectorHash> &column_sets,
                                         VectorMap<node_index, size_t> &node_to_cols,
                                         std::vector<CoordinateSet> &label_coords) const {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(
        canonical_ ? &canonical_->get_graph() : &graph_
    );

    VectorSet<Row> queued_rows;
    std::vector<node_index> queued_nodes;
    tsl::hopscotch_set<node_index> dummy_nodes;

    std::vector<std::pair<node_index, node_index>> add_base_annot;

    std::function<void(node_index, node_index)> queue_node
        = [](node_index, node_index) {};
    if (canonical_) {
        queue_node = [&](node_index node, node_index base_node) {
            assert(base_node);
            auto find_base = node_to_cols.find(base_node);
            auto row = AnnotatedDBG::graph_to_anno_index(base_node);
            if (find_base != node_to_cols.end()) {
                assert(find_base->second != nannot || queued_rows.count(row));
                return;
            }

            if (dbg_succ && !dbg_succ->get_mask()
                    && dbg_succ->get_boss().is_dummy(base_node)) {
                assert(!node_to_cols.count(node));
                dummy_nodes.emplace(node);
                return;
            }

            if (queued_rows.emplace(row).second)
                node_to_cols.emplace(base_node, nannot);
        };
    } else if (graph_.get_mode() == DeBruijnGraph::BASIC) {
        queue_node = [&](node_index node, node_index) {
            assert(node);
            auto find = node_to_cols.find(node);
            auto row = AnnotatedDBG::graph_to_anno_index(node);
            if (find != node_to_cols.end()) {
                assert(find->second != nannot || queued_rows.count(row));
                return;
            }

            if (dbg_succ && !dbg_succ->get_mask()
                    && dbg_succ->get_boss().is_dummy(node)) {
                assert(!node_to_cols.count(node));
                dummy_nodes.emplace(node);
                return;
            }

            if (queued_rows.emplace(row).second)
                node_to_cols.emplace(node, nannot);
        };
    } else {
        assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);
        queue_node = [&](node_index node, node_index base_node) {
            assert(node);
            if (base_node) {
                auto find_base = node_to_cols.find(base_node);
                if (find_base == node_to_cols.end()) {
                    auto row = AnnotatedDBG::graph_to_anno_index(base_node);
                    if (queued_rows.emplace(row).second) {
                        node_to_cols.emplace(base_node, nannot);
                        node_to_cols.emplace(node, nannot);
                        queued_nodes.emplace_back(node);
                    }
                } else if (node != base_node) {
                    node_to_cols.try_emplace(node, find_base->second);
                    if (find_base->second == nannot)
                        add_base_annot.emplace_back(node, base_node);
                }
            } else {
                assert(dbg_succ);
                assert(!dbg_succ->get_mask());
                assert(dbg_succ->get_boss().is_dummy(node));
                assert(!node_to_cols.count(node));
                dummy_nodes.emplace(node);
            }
        };
    }

    for (const auto &path : queued_paths_) {
        if (canonical_) {
            for (node_index node : path) {
                queue_node(node, canonical_->get_base_node(node));
            }
        } else if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            for (node_index node : path) {
                queue_node(node, node);
            }
        } else {
            // TODO: avoid this spelling
            std::string spelling = spell_path(graph_, path);
            auto it = path.begin();
            for (node_index base_node : map_to_nodes(graph_, spelling)) {
                assert(it != path.end());
                queue_node(*it, base_node);
                ++it;
            }
            assert(it == path.end());
        }
    }

    tsl::hopscotch_set<node_index> annotated_nodes;
    tsl::hopscotch_map<node_index, std::vector<node_index>> parents;
    for (node_index node : dummy_nodes) {
        assert(dbg_succ);
        assert(!dbg_succ->get_mask());

        // if we already discovered this via another node, move on
        node_index base_node = canonical_ ? canonical_->get_base_node(node) : node;
        assert(base_node);
        if (node_to_cols.count(base_node))
            continue;

        std::vector<std::pair<node_index, size_t>> traversal;
        std::string spelling = graph_.get_node_sequence(node);
        assert(spelling.back() != boss::BOSS::kSentinel);
        traversal.emplace_back(node, spelling.find_last_of(boss::BOSS::kSentinel) + 1);
        assert(traversal.back().second < spelling.size());

        while (traversal.size()) {
            node_index cur_node = traversal.back().first;
            size_t num_sentinels_left = traversal.back().second;
            traversal.pop_back();

            node_index cur_base_node = canonical_
                ? canonical_->get_base_node(cur_node)
                : cur_node;
            assert(cur_base_node);

            assert(dbg_succ->kmer_to_boss_index(cur_base_node) == cur_base_node);
            assert(!num_sentinels_left
                == !dbg_succ->get_boss().is_dummy(cur_base_node));

            auto find_base = node_to_cols.find(cur_base_node);
            if (find_base != node_to_cols.end()) {
                assert(canonical_ || node_to_cols.count(cur_node));

                if (!num_sentinels_left) {
                    assert(find_base->second != nannot
                        || queued_rows.count(AnnotatedDBG::graph_to_anno_index(cur_base_node)));

                    annotated_nodes.emplace(cur_node);
                }

                continue;
            }

            if (!num_sentinels_left) {
                queue_node(cur_node, cur_base_node);
                assert(node_to_cols.count(cur_base_node));
                annotated_nodes.emplace(cur_node);
                continue;
            }

            node_to_cols.try_emplace(cur_base_node, nannot);
            if (!canonical_ && cur_node != cur_base_node) {
                assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);
                node_to_cols.try_emplace(cur_node, nannot);
            }

            --num_sentinels_left;
            graph_.adjacent_outgoing_nodes(cur_node, [&](node_index next) {
                assert(graph_.get_node_sequence(next).back() != boss::BOSS::kSentinel);
                parents[next].emplace_back(cur_node);
                traversal.emplace_back(next, num_sentinels_left);
            });
        }
    }

    dummy_nodes.clear();

    auto push_node_labels = [&](node_index node,
                                auto row,
                                auto&& labels,
                                const CoordinateSet coords = {}) {
        auto do_push = [&](auto find, size_t labels_i) {
            find.value() = labels_i;
            if (has_coordinates()) {
                assert(coords.size());
                size_t coord_idx = find - node_to_cols.begin();
                if (coord_idx == label_coords.size()) {
                    label_coords.emplace_back(coords);
                } else {
                    label_coords.resize(std::max(label_coords.size(), coord_idx + 1));
                    label_coords[coord_idx] = coords;
                }
            }
        };

        node_index base_node = AnnotatedDBG::anno_to_graph_index(row);
        auto find_base = node_to_cols.find(base_node);
        assert(find_base != node_to_cols.end());
        size_t labels_i = cache_column_set_impl(column_sets, std::move(labels));
        do_push(find_base, labels_i);

        if (canonical_ || graph_.get_mode() == DeBruijnGraph::BASIC)
            return;

        if (node != base_node) {
            auto find = node_to_cols.find(node);
            assert(find != node_to_cols.end());
            assert(find->second == nannot);
            assert(find_base->second != nannot);
            do_push(find, labels_i);
        }

        if (!canonical_ && graph_.get_mode() == DeBruijnGraph::CANONICAL && base_node == node) {
            auto spelling = graph_.get_node_sequence(node);
            reverse_complement(spelling.begin(), spelling.end());
            if (node_index rc_node = map_to_nodes_sequentially(graph_, spelling)[0])
                do_push(node_to_cols.try_emplace(rc_node, nannot).first, labels_i);
        }
    };

    auto row_it = queued_rows.begin();
    auto node_it = queued_nodes.begin();
    if (has_coordinates()) {
        assert(multi_int_);
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows.values_container())) {
            assert(row_it != queued_rows.end());
            assert(!dbg_succ || dbg_succ->get_mask()
                || !dbg_succ->get_boss().is_dummy(AnnotatedDBG::anno_to_graph_index(*row_it)));
            assert(row_tuples.size());
            std::sort(row_tuples.begin(), row_tuples.end(), utils::LessFirst());
            Columns labels;
            CoordinateSet coords;
            labels.reserve(row_tuples.size());
            coords.reserve(row_tuples.size());
            for (auto&& [label, cur_coords] : row_tuples) {
                labels.push_back(label);
                coords.emplace_back(cur_coords.begin(), cur_coords.end());
            }

            assert(row_it != queued_rows.end());
            if (queued_nodes.size()) {
                assert(node_it != queued_nodes.end());
                push_node_labels(*node_it, *row_it, std::move(labels), coords);
                ++node_it;
            } else {
                push_node_labels(AnnotatedDBG::anno_to_graph_index(*row_it),
                                 *row_it, std::move(labels), coords);
            }
            ++row_it;
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows.values_container())) {
            assert(row_it != queued_rows.end());
            assert(!dbg_succ || dbg_succ->get_mask()
                || !dbg_succ->get_boss().is_dummy(AnnotatedDBG::anno_to_graph_index(*row_it)));
            if (labels.empty()) {
                logger->error("Failed\t{}:{}", AnnotatedDBG::anno_to_graph_index(*row_it),graph_.get_node_sequence(AnnotatedDBG::anno_to_graph_index(*row_it)));
            }
            assert(labels.size());
            std::sort(labels.begin(), labels.end());
            if (queued_nodes.size()) {
                assert(!canonical_ && graph_.get_mode() == DeBruijnGraph::CANONICAL);
                assert(node_it != queued_nodes.end());
                push_node_labels(*node_it, *row_it, std::move(labels));
                ++node_it;
            } else {
                push_node_labels(AnnotatedDBG::anno_to_graph_index(*row_it),
                                 *row_it, std::move(labels));
            }
            ++row_it;
        }
    }

    assert(row_it == queued_rows.end());
    assert(node_it == queued_nodes.end());

    for (const auto &[node, base_node] : add_base_annot) {
        auto find_base = node_to_cols.find(base_node);
        assert(find_base != node_to_cols.end());
        assert(find_base->second != nannot);

        auto find = node_to_cols.find(node);
        assert(find != node_to_cols.end());
        assert(find->second == nannot || find->second == find_base->second);
        find.value() = find_base->second;
        if (has_coordinates()) {
            size_t base_coord_idx = find_base - node_to_cols.begin();
            assert(base_coord_idx < label_coords.size());

            const auto &coords = label_coords[base_coord_idx];

            size_t coord_idx = find - node_to_cols.begin();
            if (coord_idx == label_coords.size()) {
                label_coords.emplace_back(coords);
            } else {
                label_coords.resize(std::max(label_coords.size(), coord_idx + 1));
                label_coords[coord_idx] = coords;
            }
        }
    }

    for (node_index node : annotated_nodes) {
        assert(parents.count(node));
        assert(get_labels(node));
        std::vector<node_index> back_traversal;
        back_traversal.emplace_back(node);
        while (back_traversal.size()) {
            node_index node = back_traversal.back();
            back_traversal.pop_back();
            assert(parents.count(node));

            auto [labels, coords] = get_labels_and_coords(node);
            assert(labels);
            assert(labels->size());

            for (node_index prev : parents[node]) {
                node_index base_node = canonical_ ? canonical_->get_base_node(prev) : prev;
                assert(canonical_ || node_to_cols.count(prev));
                assert(node_to_cols.count(base_node));
                auto [prev_labels, prev_coords] = get_labels_and_coords(prev);
                CoordinateSet merged_prev_coords;
                if (!prev_labels) {
                    if (has_coordinates()) {
                        assert(coords);
                        merged_prev_coords.reserve(coords->size());
                        for (auto &tuple : *coords) {
                            auto &prev_tuple = merged_prev_coords.emplace_back();
                            prev_tuple.reserve(tuple.size());
                            for (auto c : tuple) {
                                prev_tuple.emplace_back(c - 1);
                            }
                        }
                    }

                    push_node_labels(prev,
                                     AnnotatedDBG::graph_to_anno_index(base_node),
                                     decltype(*labels)(*labels),
                                     merged_prev_coords);
                } else {
                    Columns merged_columns;
                    if (has_coordinates()) {
                        assert(coords);
                        assert(prev_coords);
                        utils::match_indexed_values(labels->begin(), labels->end(),
                                                    coords->begin(),
                                                    prev_labels->begin(), prev_labels->end(),
                                                    prev_coords->begin(),
                                                    [&](const auto label,
                                                        const auto &c1,
                                                        const auto &c2) {
                            merged_columns.emplace_back(label);
                            auto &merge_coords = merged_prev_coords.emplace_back();
                            utils::set_union(c2.begin(), c2.end(), c1.begin(), c1.end(),
                                             std::back_inserter(merge_coords), -1);
                        });
                    } else {
                        std::set_union(labels->begin(), labels->end(),
                                       prev_labels->begin(), prev_labels->end(),
                                       std::back_inserter(merged_columns));
                    }

                    push_node_labels(prev,
                                     AnnotatedDBG::graph_to_anno_index(base_node),
                                     std::move(merged_columns),
                                     merged_prev_coords);
                }

                if (parents.count(prev)) {
                    assert(get_labels(prev));
                    back_traversal.emplace_back(prev);
                }
            }
        }
    }

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols) {
        assert(val != nannot);
    }
#endif

}

auto AnnotationBuffer::get_labels_and_coords(node_index node) const
        -> std::pair<const Columns*, const CoordinateSet*> {
    std::pair<const Columns*, const CoordinateSet*> ret_val { nullptr, nullptr };
    if (!node) {
        ret_val.first = &column_sets_.data()[0];
        return ret_val;
    }

    if (canonical_)
        node = canonical_->get_base_node(node);

    auto it = node_to_cols_.find(node);

    // if the node hasn't been seen before, or if its annotations haven't
    // been fetched, return nothing
    if (it == node_to_cols_.end() || it->second == nannot)
        return ret_val;

    ret_val.first = &column_sets_.data()[it->second];

    if (has_coordinates()) {
        assert(static_cast<size_t>(it - node_to_cols_.begin()) < label_coords_.size());
        ret_val.second = &label_coords_[it - node_to_cols_.begin()];
    }

    return ret_val;
}

bool AnnotationBuffer::allow_label_change() const {
    return graph_.get_extension_threadsafe<HLLWrapper<>>();
}

Alignment::score_t AnnotationBuffer
::get_label_change_score(Alignment::Column col_a, Alignment::Column col_b) const {
    if (col_a == col_b)
        return 0;

    Alignment::score_t label_change_score = DBGAlignerConfig::ninf;

    if (const auto *hll_wrapper = graph_.get_extension_threadsafe<HLLWrapper<>>()) {
        const auto &hll = hll_wrapper->data();
        uint64_t a_size = hll.num_relations_in_column(col_a);
        uint64_t b_size = hll.num_relations_in_column(col_b);

        double dbsize = b_size;
        double logdbsize = log2(dbsize);

        auto &cache_a = cache_[col_a];
        auto [cache_find, cache_inserted] = cache_a.try_emplace(col_b, 0.0);

        double union_size;
        if (cache_inserted) {
            union_size = hll.estimate_column_union_cardinality(col_a, col_b);
            cache_find.value() = union_size;
        } else {
            union_size = cache_find->second;
        }

        uint64_t size_sum = a_size + b_size;
        if (union_size < size_sum)
            label_change_score = log2(std::min(dbsize, size_sum - union_size)) - logdbsize;
    }

    assert(label_change_score <= 0);
    return label_change_score;
}

} // namespace align
} // namespace graph
} // namespace mtg
