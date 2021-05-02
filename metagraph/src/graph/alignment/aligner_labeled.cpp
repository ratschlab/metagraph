#include "aligner_labeled.hpp"

#include <tsl/hopscotch_set.h>

#include "aligner_prefix_suffix.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/utils/template_utils.hpp"

namespace mtg {
namespace graph {
namespace align {


bool check_targets(const AnnotatedDBG &anno_graph,
                   const Alignment<DeBruijnGraph::node_index> &path) {
    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<uint64_t> ref_targets;
    for (const std::string &label : anno_graph.get_labels(path.get_sequence(), 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (uint64_t target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}


void
process_seq_path(const DeBruijnGraph &graph,
                 std::string_view query,
                 const std::vector<DeBruijnGraph::node_index> &query_nodes,
                 const std::function<void(AnnotatedDBG::row_index, size_t)> &callback) {
    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    if (canonical) {
        if (query_nodes.size()) {
            auto first = std::find_if(query_nodes.begin(), query_nodes.end(),
                                      [](auto i) -> bool { return i; });
            if (first == query_nodes.end())
                return;

            size_t start = first - query_nodes.begin();

            if (canonical->get_base_node(*first) == *first) {
                for (size_t i = start; i < query_nodes.size(); ++i) {
                    if (query_nodes[i] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i])
                        ), i);
                    }
                }
            } else {
                for (size_t i = query_nodes.size(); i > start; --i) {
                    if (query_nodes[i - 1] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i - 1])
                        ), i - 1);
                    }
                }
            }
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos)
                callback(AnnotatedDBG::graph_to_anno_index(query_nodes[i]), i);
        }
    } else {
        size_t i = 0;
        if (query.front() == '#') {
            std::string map_query
                = graph.get_node_sequence(query_nodes[0]).substr(0, graph.get_k());
            map_query += query.substr(graph.get_k());
            graph.map_to_nodes(map_query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        } else {
            graph.map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        }
        assert(i == query_nodes.size());
    }
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>::init_backtrack() const {
    std::vector<uint64_t> added_rows;
    std::vector<node_index> added_nodes;

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&this->graph_);
    if (!dbg_succ) {
        if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_))
            dbg_succ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
    }

    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    DefaultColumnExtender<NodeType>::call_visited_nodes([&](node_index node, size_t, size_t) {
        auto [it, inserted] = targets_.emplace(node, 0);
        if (inserted) {
            process_seq_path(this->graph_, std::string(this->graph_.get_k(), '#'),
                             { node }, [&](auto row, size_t) {
                if (!boss || boss->get_W(dbg_succ->kmer_to_boss_index(AnnotatedDBG::anno_to_graph_index(row)))) {
                    added_rows.push_back(row);
                    added_nodes.push_back(node);
                }
            });
        }
    });

    auto it = added_nodes.begin();
    for (auto &labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows)) {
        assert(it != added_nodes.end());
        auto jt = targets_set_.emplace(labels).first;
        targets_[*it] = jt - targets_set_.begin();
        ++it;
    }
    assert(it == added_nodes.end());
}

template <typename NodeType>
auto LabeledBacktrackingExtender<NodeType>
::backtrack(score_t min_path_score,
            AlignNode best_node,
            tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
            std::vector<DBGAlignment> &extensions) const -> std::vector<AlignNode> {
    size_t target_id = targets_[std::get<0>(best_node)];
    if (!target_id)
        return {};

    std::vector<DBGAlignment> next_extension;
    tsl::hopscotch_set<AlignNode, AlignNodeHash> dummy;
    auto track = DefaultColumnExtender<NodeType>::backtrack(min_path_score, best_node,
                                                            dummy, next_extension);

    assert(next_extension.size() <= 1);
    if (next_extension.empty() || next_extension[0].get_offset()) {
        prev_starts.insert(dummy.begin(), dummy.end());
        return track;
    }

    dummy = tsl::hopscotch_set<AlignNode, AlignNodeHash>();

    DBGAlignment alignment = std::move(next_extension[0]);
    assert(alignment.is_valid(this->graph_, &this->config_));
    assert(alignment.back() == std::get<0>(best_node));

    Vector<uint64_t> target_intersection = *(targets_set_.begin() + target_id);

    AlignmentSuffix<node_index> suffix(alignment, this->config_, this->graph_);
    while (!suffix.get_added_offset())
        ++suffix;

    auto suffix_shift = [&suffix]() {
        assert(!suffix.reof());
        --suffix;
        while (!suffix.reof() && suffix.get_front_op() == Cigar::INSERTION) {
            --suffix;
        }
    };

    suffix_shift();

    auto it = track.begin();
    assert(it != track.end());

    prev_starts.emplace(*it);
    ++it;

    for (size_t i = alignment.size() - 1; i > 0; --i) {
        const auto &cur_targets = *(targets_set_.begin() + targets_[alignment[i - 1]]);
        Vector<uint64_t> inter;
        std::set_intersection(target_intersection.begin(), target_intersection.end(),
                              cur_targets.begin(), cur_targets.end(),
                              std::back_inserter(inter));

        if (inter.empty())
            break;

        if (it != track.end()) {
            Vector<uint64_t> diff;
            std::set_difference(cur_targets.begin(), cur_targets.end(),
                                target_intersection.begin(), target_intersection.end(),
                                std::back_inserter(diff));
            if (inter.size() == target_intersection.size() || diff.empty()) {
                prev_starts.emplace(*it);
                ++it;
            } else {
                it = track.end();
            }
        }

        if (inter.size() < target_intersection.size()) {
            if (suffix.get_front_op() != Cigar::DELETION) {
                extensions.emplace_back(suffix);
                extensions.back().target_columns = target_intersection;
                assert(check_targets(anno_graph_, extensions.back()));
            }
        } else {
            --suffix;
            if (!suffix.reof() && suffix.get_front_op() == Cigar::DELETION) {
                ++suffix;
                if (suffix.get_front_op() != Cigar::DELETION) {
                    extensions.emplace_back(suffix);
                    extensions.back().target_columns = target_intersection;
                    assert(check_targets(anno_graph_, extensions.back()));
                }
            } else {
                ++suffix;
            }
        }

        suffix_shift();
        std::swap(target_intersection, inter);
    }

    if (target_intersection.size() && (suffix.reof() || suffix.get_front_op() != Cigar::DELETION)) {
        extensions.emplace_back(suffix);
        extensions.back().target_columns = target_intersection;
        assert(check_targets(anno_graph_, extensions.back()));
    }

    return track;
}

template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
