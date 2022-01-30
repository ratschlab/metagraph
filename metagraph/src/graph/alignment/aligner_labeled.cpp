#include "aligner_labeled.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_lcs.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

using MIM = annot::matrix::MultiIntMatrix;


AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const MIM*>(&annotator_.get_matrix())),
        labels_set_({ {} }) {
    if (multi_int_ && graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        common::logger->warn("Coordinates not supported when aligning to CANONICAL "
                             "or PRIMARY mode graphs");
    }
}

void AnnotationBuffer::flush() {
    if (added_rows_.empty())
        return;

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels) {
        assert(node_it != added_nodes_.end());
        assert(row_it != added_rows_.end());
        auto label_it = labels_set_.emplace(std::forward<decltype(labels)>(labels)).first;
        assert(labels_.count(*node_it));
        assert(labels_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));
        assert(labels_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        labels_[*node_it].second = label_it - labels_set_.begin();
        labels_[AnnotatedDBG::anno_to_graph_index(*row_it)].second = label_it - labels_set_.begin();
    };

    auto node_it = added_nodes_.begin();
    auto row_it = added_rows_.begin();
    if (multi_int_) {
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(added_rows_)) {
            Vector<Column> labels;
            labels.reserve(row_tuples.size());
            label_coords_.emplace_back();
            label_coords_.back().reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                label_coords_.back().emplace_back(std::forward<decltype(coords)>(coords));
            }

            push_node_labels(node_it++, row_it++, std::move(labels));
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(added_rows_)) {
            push_node_labels(node_it++, row_it++, std::forward<decltype(labels)>(labels));
        }
    }

    assert(node_it == added_nodes_.end());
    assert(row_it == added_rows_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

struct HasNext {
    HasNext(ssize_t diff) : diff_(diff) {}

    template <class InIt1, class InIt2>
    bool operator()(InIt1 a_begin, InIt1 a_end, InIt2 b_begin, InIt2 b_end) const {
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + diff_ < *b_begin) {
                ++a_begin;
            } else if (*a_begin + diff_ > *b_begin) {
                ++b_begin;
            } else {
                return true;
            }
        }

        return false;
    }

    ssize_t diff_;
};

template <class AIt, class BIt, class OutIt, class OutIt2>
void set_intersection_difference(AIt a_begin,
                                 AIt a_end,
                                 BIt b_begin,
                                 BIt b_end,
                                 OutIt intersection_out,
                                 OutIt2 diff_out) {
    while (a_begin != a_end) {
        if (b_begin == b_end || *a_begin < *b_begin) {
            *diff_out = *a_begin;
            ++diff_out;
            ++a_begin;
        } else if (*a_begin > *b_begin) {
            ++b_begin;
        } else {
            *intersection_out = *a_begin;
            ++intersection_out;
            ++a_begin;
            ++b_begin;
        }
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (DefaultColumnExtender::set_seed(seed)) {
        node_labels_.resize(1);
        label_changed_.resize(1);
        if (seed.label_columns.size()) {
            if (auto labels = labeled_graph_.get_labels_and_coordinates(seed.get_nodes()[0]).first) {
                Vector<Column> inter;
                std::set_intersection(seed.label_columns.begin(),
                                      seed.label_columns.end(),
                                      labels->get().begin(), labels->get().end(),
                                      std::back_inserter(inter));
                node_labels_[0] = std::move(inter);
            } else {
                node_labels_[0] = seed.label_columns;
            }
        } else {
            node_labels_[0] = std::nullopt;
        }

        return true;
    }

    return false;
}

void LabeledExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char /* last char */, score_t)> &callback,
                size_t table_i,
                bool /* force_fixed_seed */) {
    assert(label_changed_.size() == node_labels_.size());
    assert(node_labels_.size() == table.size());
    size_t next_offset = std::get<6>(table[table_i]) + 1;
    bool in_seed = next_offset - this->seed_->get_offset()
                    < this->seed_->get_sequence().size();
    assert(node == std::get<3>(table[table_i]));

    if (in_seed && next_offset < graph_->get_k()) {
        DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                             [&](node_index next, char c, score_t score) {
            callback(next, c, score);
            node_labels_.emplace_back(node_labels_[0]);
            label_changed_.emplace_back(false);
        }, table_i);

        return;
    }

    assert(graph_->get_node_sequence(node).back() == std::get<5>(table[table_i]));
    for ( ; last_buffered_table_i_ < table.size(); ++last_buffered_table_i_) {
        labeled_graph_.add_node(std::get<3>(table[last_buffered_table_i_]));
    }

    std::vector<std::tuple<node_index, char, score_t, node_index>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                         [&](node_index next, char c, score_t score) {
        outgoing.emplace_back(next, c, score, labeled_graph_.add_node(next));
    }, table_i);

    // if (outgoing.size() > 1)
    //     labeled_graph_.flush();

    // const std::optional<Vector<Column>> &base_labels = outgoing.size() > 1
    //     ? node_labels_[table_i]
    //     : node_labels_[std::get<10>(table[table_i])];
    labeled_graph_.flush();
    const std::optional<Vector<Column>> &base_labels = node_labels_[table_i];

    const DBGSuccinct *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph_->get_base_graph());

    size_t expand_seed = dbg_succ
        ? dbg_succ->get_boss().get_k()
            - std::min(dbg_succ->get_boss().get_k(), config_.label_change_search_width)
        : graph_->get_k();

    std::shared_ptr<NodeLCS> lcs = expand_seed < graph_->get_k()
        ? dbg_succ->get_extension<NodeLCS>()
        : nullptr;

    std::vector<std::tuple<Vector<Column>, Vector<Column>, double>> inter_diff(outgoing.size());
    size_t found_diff = 0;
    size_t found_subset = 0;
    for (size_t i = 0; i < outgoing.size(); ++i) {
        const auto &[next, c, score, next_base] = outgoing[i];
        if (auto next_labels = labeled_graph_.get_labels_and_coordinates(next).first) {
            auto &[inter, diff, lclogprob] = inter_diff[i];
            if (base_labels) {
                set_intersection_difference(next_labels->get().begin(),
                                            next_labels->get().end(),
                                            base_labels->begin(), base_labels->end(),
                                            std::back_inserter(inter),
                                            std::back_inserter(diff));

                if (inter.size() || !lcs) {
                    ++found_subset;
                } else {
                    ++found_diff;
                }

            } else {
                inter = next_labels->get();
                ++found_subset;
            }
        } else {
            ++found_subset;
        }
    }

    assert(found_diff + found_subset == outgoing.size());
    double sum_diff_probs = 0.0;

    if (found_diff) {
        assert(dbg_succ);
        assert(lcs);

        const boss::BOSS &boss = dbg_succ->get_boss();

        node_index node_base = labeled_graph_.get_base_node(node);
        assert(node_base);
        boss::BOSS::edge_index edge_base = dbg_succ->kmer_to_boss_index(node_base);
        boss::BOSS::edge_index last_edge_base = boss.succ_last(edge_base);
        std::pair<boss::BOSS::edge_index, boss::BOSS::edge_index> node_range {
            boss.pred_last(last_edge_base - 1) + 1, last_edge_base
        };
        node_range = lcs->expand(node_range.first, node_range.second,
                                 boss.get_k(), boss.get_k() - expand_seed);
        assert((*lcs)[node_range.first] < expand_seed);
        assert(node_range.second + 1 == lcs->data().size()
            || (*lcs)[node_range.second + 1] < expand_seed);

        std::pair<boss::BOSS::edge_index, boss::BOSS::edge_index> next_range;
        if (node == node_base && graph_->is_base_node(node_base)) {
            // $$$ATGC X -> $$ATGCX
            next_range = node_range;
            assert(boss.encode(std::get<5>(table[table_i]))
                    == boss.get_W(dbg_succ->kmer_to_boss_index(node)) % boss.alph_size);
            bool check = boss.tighten_range(&next_range.first, &next_range.second,
                                            boss.encode(std::get<5>(table[table_i])));
            std::ignore = check;
            assert(check);
        } else {
            // $$$ATGC X -> $$ZYATG C
            next_range.second = boss.succ_last(boss.bwd(last_edge_base));
            next_range.first = boss.pred_last(next_range.second - 1) + 1;
            if (expand_seed + 2 < boss.get_k()) {
                next_range = lcs->expand(next_range.first, next_range.second,
                                         boss.get_k(), boss.get_k() - (expand_seed + 2));
            }
        }

        for (auto i = node_range.first; i <= node_range.second; ++i) {
            if (node_index n = dbg_succ->boss_to_kmer_index(i))
                labeled_graph_.add_node(n);
        }

        for (auto i = next_range.first; i <= next_range.second; ++i) {
            if (node_index n = dbg_succ->boss_to_kmer_index(i))
                labeled_graph_.add_node(n);
        }

        labeled_graph_.flush();

        boss::BOSS::TAlphabet cur_edge_label = boss.get_W(edge_base) % boss.alph_size;

        for (size_t j = 0; j < inter_diff.size(); ++j) {
            auto &[inter, diff, lclogprob] = inter_diff[j];
            if (inter.size() || diff.empty() || !found_diff)
                continue;

#ifndef NDEBUG
            const auto &[next, c, score, next_base] = outgoing[j];
            assert(next != next_base || node != node_base || (
                next_range.first <= dbg_succ->kmer_to_boss_index(next)
                && next_range.second >= dbg_succ->kmer_to_boss_index(next)
            ));
#endif

            /**
             * IL: independent label
             * score = log2(P(traversing IL)) - log2(P(IL in node) * P(IL in next))
             *       = log2(P(traversing IL | IL in node, IL in next))
             */
            tsl::hopscotch_map<Column, std::pair<size_t, size_t>> shared_label_counts;
            tsl::hopscotch_set<Column> exclude { inter.begin(), inter.end() };
            for (Column c : diff) {
                exclude.emplace(c);
            }

            size_t node_inc_count = 0;
            size_t node_exc_count = 0;
            auto update_node_count = [&](boss::BOSS::edge_index i) {
                if (node_index n = dbg_succ->boss_to_kmer_index(i)) {
                    if (auto labels = labeled_graph_.get_labels_and_coordinates(n).first) {
                        for (Column col : labels->get()) {
                            if (!exclude.count(col)) {
                                ++shared_label_counts[col].first;
                                ++node_exc_count;
                            } else {
                                ++node_inc_count;
                            }
                        }
                    }
                }
            };
            for (auto i = boss.succ_W(node_range.first, cur_edge_label);
                    i <= node_range.second;
                    i = boss.succ_W(i + 1, cur_edge_label)) {
                update_node_count(i);
            }
            for (auto i = boss.succ_W(node_range.first, cur_edge_label + boss.alph_size);
                    i <= node_range.second;
                    i = boss.succ_W(i + 1, cur_edge_label + boss.alph_size)) {
                update_node_count(i);
            }

            size_t next_inc_count = 0;
            size_t next_exc_count = 0;
            for (auto i = next_range.first; i <= next_range.second; ++i) {
                if (node_index n = dbg_succ->boss_to_kmer_index(i)) {
                    if (auto labels = labeled_graph_.get_labels_and_coordinates(n).first) {
                        for (Column col : labels->get()) {
                            if (!exclude.count(col)) {
                                ++shared_label_counts[col].second;
                                ++next_exc_count;
                            } else {
                                ++next_inc_count;
                            }
                        }
                    }
                }
            }

            if (shared_label_counts.empty()) {
                lclogprob = config_.ninf;
                DEBUG_LOG("Label change score: {}", lclogprob);
                continue;
            }

            constexpr double pseudocount = 0.1;
            double numer = pseudocount * pseudocount * shared_label_counts.size();
            for (const auto &[col, count_pair] : shared_label_counts) {
                assert(count_pair.first || count_pair.second);
                numer += pseudocount * (count_pair.first + count_pair.second)
                            + count_pair.first * count_pair.second;
            }
            assert(numer > 0);

            double denom_pseudocount = pseudocount * shared_label_counts.size();
            double node_total_count = denom_pseudocount + node_inc_count + node_exc_count;
            double next_total_count = denom_pseudocount + next_inc_count + next_exc_count;
            lclogprob = log2(numer) - log2(node_total_count) - log2(next_total_count);
            assert(lclogprob < 0);

            sum_diff_probs += 1.0 - std::pow(2.0, lclogprob);
            // std::cerr << "test\t" << numer << "\t" << node_total_count << "\t" << next_total_count << "\t" << lclogprob << "\t" << outgoing.size();
            lclogprob = std::floor(lclogprob - log2(static_cast<double>(outgoing.size())));
            // std::cerr << "\tlabel c\t" << lclogprob << "\n";
            DEBUG_LOG("Label change score: {}", lclogprob);
        }

        sum_diff_probs = std::floor(log2((sum_diff_probs + 1.0)
                            / static_cast<double>(outgoing.size())));
        // std::cerr << "\tlabel k\t" << sum_diff_probs << "\n";
        DEBUG_LOG("Label preservation score: {}", sum_diff_probs);
    }

    for (size_t i = 0; i < inter_diff.size(); ++i) {
        const auto &[next, c, score, next_base] = outgoing[i];
        auto &[inter, diff, lclogprob] = inter_diff[i];
        if (diff.empty()) {
            label_changed_.emplace_back(false);
            node_labels_.emplace_back(std::move(inter));
            callback(next, c, score + sum_diff_probs);
        } else if (!found_diff || lclogprob > config_.ninf) {
            label_changed_.emplace_back(true);
            node_labels_.emplace_back(std::move(diff));
            callback(next, c, score + (found_diff ? lclogprob : -config_.extra_penalty));
        }
    }

    // assert(node_labels_.size() == table.size() + outgoing.size());
}

void LabeledExtender
::call_alignments(score_t cur_cell_score,
                  score_t end_score,
                  score_t min_path_score,
                  const std::vector<node_index> &path,
                  const std::vector<size_t> &trace,
                  size_t table_i,
                  const Cigar &ops,
                  size_t clipping,
                  size_t offset,
                  std::string_view window,
                  const std::string &match,
                  score_t extra_penalty,
                  const std::function<void(Alignment&&)> &callback) {
    DefaultColumnExtender::call_alignments(cur_cell_score, end_score, min_path_score,
                                           path, trace, table_i, ops, clipping,
                                           offset, window, match, extra_penalty,
                                           [&](Alignment&& aln) {
        // store the label and coordinate information
        aln.label_encoder = &labeled_graph_.get_annotator().get_label_encoder();

        auto get_annotation = [&](size_t table_i, node_index node) -> Vector<Column> {
            if (std::optional<Vector<Column>> fetch_labels = node_labels_[table_i])
                return *fetch_labels;

            auto labels = labeled_graph_.get_labels_and_coordinates(node).first;
            assert(labels);
            size_t last_fork_i = std::get<10>(this->table[table_i]);
            if (std::optional<Vector<Column>> fetch_labels = node_labels_[last_fork_i]) {
                Vector<Column> inter;
                std::set_intersection(labels->get().begin(), labels->get().end(),
                                      fetch_labels->begin(), fetch_labels->end(),
                                      std::back_inserter(inter));
                return inter;
            }

            return labels->get();
        };

        aln.label_columns = get_annotation(trace[0], aln.get_nodes().back());

        if (extra_penalty) {
            auto it = trace.begin() + 1;
            for (auto jt = aln.get_nodes().rbegin() + 1;
                    jt != aln.get_nodes().rend(); ++it, ++jt) {
                if (label_changed_[*(it - 1)]) {
                    auto next = get_annotation(*it, *jt);
                    Vector<Column> label_union;
                    std::set_union(aln.label_columns.begin(), aln.label_columns.end(),
                                   next.begin(), next.end(), std::back_inserter(label_union));
                    std::swap(label_union, aln.label_columns);
                }
            }
        }

        callback(std::move(aln));
    });
}

void LabeledBacktrackingExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char /* last char */, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    bool in_seed = std::get<6>(table[table_i]) + 1 - this->seed_->get_offset()
                        < this->seed_->get_sequence().size();
    size_t next_offset = std::get<6>(table[table_i]) + 1;
    assert(node == std::get<3>(table[table_i]));

    if (in_seed && (next_offset < graph_->get_k() || force_fixed_seed || this->fixed_seed())) {
        DefaultColumnExtender::call_outgoing(node, max_prefetch_distance, callback,
                                             table_i, force_fixed_seed);
        return;
    }

    assert(graph_->get_node_sequence(node).back() == std::get<5>(table[table_i]));
    for ( ; last_buffered_table_i_ < table.size(); ++last_buffered_table_i_) {
        labeled_graph_.add_node(std::get<3>(table[last_buffered_table_i_]));
    }

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                         [&](node_index next, char c, score_t score) {
        outgoing.emplace_back(next, c, score);
        labeled_graph_.add_node(next);
    }, table_i, force_fixed_seed);

    if (outgoing.size() > 1)
        labeled_graph_.flush();

    const Alignment &seed = *this->seed_;

    auto [base_labels, base_coords] = labeled_graph_.get_labels_and_coordinates(node);

    if (base_coords) {
        // check label and coordinate consistency

        // first, determine a base node from which to compare coordinates
        // by default, node is used (the parent node of next)
        ssize_t dist = 1;
        bool rev_align = dynamic_cast<const RCDBG*>(graph_);
        if (seed.label_coordinates.size()) {
            // if the seed has coordinates, use the seed as the base
            base_labels = std::cref(seed.label_columns);
            base_coords = std::cref(seed.label_coordinates);
            ssize_t offset = std::get<6>(table[table_i]);
            ssize_t dist_from_origin = offset - (seed.get_offset() - 1);
            dist = dist_from_origin - seed.get_offset()
                - seed.get_sequence().size()
                + seed.get_nodes().size()
                - (seed.get_sequence().size() - graph_->get_k()) * rev_align;
        } else {
            // otherwise, check the first node in the traversal
            ssize_t k = graph_->get_k();
            std::tie(base_labels, base_coords)
                = labeled_graph_.get_labels_and_coordinates(std::get<3>(table[0]));
            ssize_t base_offset = std::get<6>(table[0]);
            if (!base_coords || base_coords->get().empty()) {
                // if the first node is not annotated, then check the last fork point
                size_t base_table_i = xdrop_cutoffs_[std::get<9>(table[table_i])].first;
                std::tie(base_labels, base_coords)
                    = labeled_graph_.get_labels_and_coordinates(std::get<3>(table[base_table_i]));
                base_offset = std::get<6>(table[base_table_i]);
            }

            if (base_coords) {
                // determine the distance of each child node to the selected base node
                dist = std::get<6>(table[table_i]) - base_offset + 1;
                if (base_offset < k)
                    dist -= k - seed.get_offset();
            }

            // if the seed has stored labels, but no stored coordinates, check
            // the base node to make sure that it has intersecting labels with the seed
            if (seed.label_columns.size()
                    && !utils::share_element(base_labels->get().begin(),
                                             base_labels->get().end(),
                                             seed.label_columns.begin(),
                                             seed.label_columns.end())) {
                return;
            }
        }

        for (const auto &[next, c, score] : outgoing) {
            if (base_labels->get().empty()) {
                callback(next, c, score);
                continue;
            }

            auto [next_labels, next_coords]
                = labeled_graph_.get_labels_and_coordinates(next);

            // check coordinate consistency later if they are not cached now
            if (!next_coords) {
                callback(next, c, score);
                continue;
            }

            // if we are traversing backwards, then negate the coordinate delta
            ssize_t dist_sign = rev_align ? -1 : 1;

            if (!utils::indexed_set_find<HasNext>(base_labels->get().begin(),
                                                  base_labels->get().end(),
                                                  base_coords->get().begin(),
                                                  next_labels->get().begin(),
                                                  next_labels->get().end(),
                                                  next_coords->get().begin(),
                                                  dist * dist_sign))
                continue;

            callback(next, c, score);
        }
    } else if (base_labels) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (const auto &[next, c, score] : outgoing) {
            if (base_labels->get().empty()) {
                callback(next, c, score);
                continue;
            }

            auto next_labels = labeled_graph_.get_labels(next);

            // if not fetched, check later
            if (!next_labels) {
                callback(next, c, score);
                continue;
            }

            if ((seed.label_columns.size()
                    && !utils::share_element(seed.label_columns.begin(),
                                             seed.label_columns.end(),
                                             next_labels->get().begin(),
                                             next_labels->get().end()))
                    || !utils::share_element(base_labels->get().begin(),
                                             base_labels->get().end(),
                                             next_labels->get().begin(),
                                             next_labels->get().end()))
                continue;

            callback(next, c, score);
        }
    } else {
        for (const auto &[next, c, score] : outgoing) {
            callback(next, c, score);
        }
    }
}

auto AnnotationBuffer::add_path(const std::vector<node_index> &path, std::string query)
        -> std::pair<std::vector<node_index>, bool> {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    if (path.empty())
        return {};

    const DeBruijnGraph &base_graph = graph_.get_base_graph();
    bool base_is_canonical = base_graph.get_mode() == DeBruijnGraph::CANONICAL;

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    if (base_is_canonical && query.front() == '#')
        query = graph_.get_node_sequence(path[0]) + query.substr(graph_.get_k());

    auto call_node = [&](node_index start_node, node_index start_base_node) {
        if (start_base_node != DeBruijnGraph::npos) {
            if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(start_base_node)))
                return; // skip dummy nodes

            std::vector<std::pair<node_index, node_index>> to_add;
            to_add.emplace_back(start_node, start_base_node);

            if (const auto *row_diff = dynamic_cast<const annot::binmat::IRowDiff*>(&annotator_.get_matrix())) {
                assert(boss);
                const auto &anchor = row_diff->anchor();
                const auto &fork_succ = row_diff->fork_succ();
                auto edge = dbg_succ->kmer_to_boss_index(start_base_node);
                while (!anchor[edge] && !labels_.count(start_base_node)) {
                    edge = boss->row_diff_successor(edge, fork_succ);
                    start_base_node = dbg_succ->boss_to_kmer_index(edge);
                    to_add.emplace_back(start_base_node, start_base_node);
                }
            }

            for (auto&& [node, base_node] : to_add) {
                Row row = AnnotatedDBG::graph_to_anno_index(base_node);
                if (labels_.emplace(node, std::make_pair(row, nannot)).second
                        && (node == base_node
                            || labels_.emplace(base_node, std::make_pair(row, nannot)).second)) {
                    added_rows_.push_back(row);
                    added_nodes_.push_back(node);
                } else {
                    auto find_n = labels_.find(node);
                    auto find_b = labels_.find(base_node);
                    assert(find_n != labels_.end());
                    assert(find_b != labels_.end());
                    assert(find_n->second.first == find_b->second.first);
                    auto label_i = std::min(find_n->second.second, find_b->second.second);
                    find_n.value().second = label_i;
                    find_b.value().second = label_i;
                }
            }
        }
    };

    auto base_path = graph_.get_base_path(path, query);
    if (!base_path.second) {
        for (size_t i = 0; i < base_path.first.size(); ++i) {
            call_node(path[i], base_path.first[i]);
        }
    } else {
        auto it = path.rbegin();
        for (size_t i = 0; i < base_path.first.size(); ++i, ++it) {
            call_node(*it, base_path.first[i]);
        }
    }

    return base_path;
}

auto AnnotationBuffer::add_node(node_index node) -> node_index {
    return add_path({ node }, std::string(graph_.get_k(), '#')).first[0];
}

bool LabeledBacktrackingExtender::skip_backtrack_start(size_t i) {
    label_intersection_.clear();
    label_intersection_coords_.clear();

    if (!this->prev_starts.emplace(i).second)
        return false;

    // if backtracking hasn't been started from here yet, get its labels
    node_index node = std::get<3>(this->table[i]);

    auto label_find = diff_label_sets_.find(i);
    if (label_find != diff_label_sets_.end()) {
        // extract a subset of the labels if this node was previously traversed
        std::tie(label_intersection_, label_intersection_coords_) = label_find->second;

    } else {
        // otherwise, take the full label set
        auto [fetch_labels, fetch_coords]
            = labeled_graph_.get_labels_and_coordinates(node);
        if (fetch_labels) {
            label_intersection_ = fetch_labels->get();
            if (fetch_coords) {
                label_intersection_coords_ = fetch_coords->get();
                if (seed_label_coordinates_.size()) {
                    const Alignment &seed = *this->seed_;
                    bool rev_align = dynamic_cast<const RCDBG*>(graph_);

                    ssize_t offset = std::get<6>(this->table[i]);
                    ssize_t dist_from_origin = offset - (seed.get_offset() - 1) - 1;
                    ssize_t dist = dist_from_origin - seed.get_offset()
                                    - seed.get_sequence().size()
                                    + seed.get_nodes().size()
                                    - (seed.get_sequence().size() - graph_->get_k()) * rev_align;
                    ssize_t dist_sign = rev_align ? 1 : -1;

                    Vector<Column> labels_inter;
                    Vector<Tuple> coords_inter;
                    utils::indexed_set_op<Tuple, CoordIntersection>(
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        seed_labels_.begin(),
                        seed_labels_.end(),
                        seed_label_coordinates_.begin(),
                        std::back_inserter(labels_inter), std::back_inserter(coords_inter),
                        dist * dist_sign
                    );
                    std::swap(labels_inter, label_intersection_);
                    std::swap(coords_inter, label_intersection_coords_);
                }
            } else if (seed_labels_.size()) {
                Vector<Column> labels_inter;
                std::set_intersection(seed_labels_.begin(),
                                      seed_labels_.end(),
                                      label_intersection_.begin(),
                                      label_intersection_.end(),
                                      std::back_inserter(labels_inter));
                std::swap(labels_inter, label_intersection_);
            }
        }
    }

    // we already have the labels for the first node in the path
    last_path_size_ = 1;

    assert(!labeled_graph_.get_coordinate_matrix()
        || label_intersection_.size() == label_intersection_coords_.size());

    cur_min_path_score_ = label_intersection_.size()
        ? extensions_.get_min_path_score(label_intersection_)
        : std::numeric_limits<score_t>::max();

    // skip backtracking from this node if no labels could be determined for it
    return label_intersection_.empty();
}

struct CoordDiff {
    CoordDiff(size_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin + offset_ < *b_begin) {
                *out = *a_begin;
                ++out;
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                ++a_begin;
                ++b_begin;
            }
        }
    }

    size_t offset_;
};

void LabeledBacktrackingExtender
::call_alignments(score_t cur_cell_score,
                  score_t end_score,
                  score_t min_path_score,
                  const std::vector<node_index> &path,
                  const std::vector<size_t> &trace,
                  size_t table_i,
                  const Cigar &ops,
                  size_t clipping,
                  size_t offset,
                  std::string_view window,
                  const std::string &match,
                  score_t extra_penalty,
                  const std::function<void(Alignment&&)> & /* callback */) {
    assert(path.size());
    assert(ops.size());
    assert(label_intersection_.size());

    score_t alignment_score = end_score - (clipping ? cur_cell_score : 0);
    if (cur_min_path_score_ > alignment_score)
        return;

    // Normally, we start from the end of the alignment and reconstruct the
    // alignment backwards (shifting coordinates by -1).
    // If we are aligning backwards, then the coordinates need to be shifted by
    // 1 instead.
    ssize_t coord_step = dynamic_cast<const RCDBG*>(graph_) ? 1 : -1;

    size_t label_path_end = trace.size()
        - std::min(this->graph_->get_k(), trace.size()) + 1;
    assert(label_path_end <= path.size());
    assert(label_path_end >= last_path_size_);
    constexpr uint64_t ncoord = std::numeric_limits<uint64_t>::max();
    for ( ; last_path_size_ < label_path_end; ++last_path_size_) {
        // update current coordinates
        for (auto &coords : label_intersection_coords_) {
            for (uint64_t &c : coords) {
                if (!c && coord_step == -1) {
                    c = ncoord;
                } else {
                    c += coord_step;
                }
            }

            // if any coordinates went out of bounds, remove them
            coords.erase(std::remove_if(coords.begin(), coords.end(),
                                        [ncoord](const auto &a) { return a == ncoord; }),
                         coords.end());
        }

        size_t old_size = 0;
        size_t cur_size = 0;
        size_t new_size = 0;

        const Vector<Column> *label_set = nullptr;
        const Vector<Tuple> *label_coord = nullptr;

        auto [fetch_labels, fetch_coords]
                = labeled_graph_.get_labels_and_coordinates(path[last_path_size_]);
        if (fetch_coords) {
            // backtrack if coordinates are consistent
            const Vector<Column> &labels = fetch_labels->get();
            const Vector<Tuple> &coords = fetch_coords->get();
            Vector<Column> label_inter;
            Vector<Tuple> coord_inter;
            label_set = &labels;
            label_coord = &coords;
            std::tie(old_size, cur_size, new_size)
                = utils::indexed_set_op<Tuple, CoordIntersection>(
                    label_intersection_.begin(), label_intersection_.end(),
                    label_intersection_coords_.begin(),
                    labels.begin(), labels.end(), coords.begin(),
                    std::back_inserter(label_inter), std::back_inserter(coord_inter)
                );

            std::swap(label_intersection_, label_inter);
            std::swap(label_intersection_coords_, coord_inter);

            if (label_intersection_.empty()) {
                cur_min_path_score_ = std::numeric_limits<score_t>::max();
                return;
            }

        } else if (fetch_labels) {
            // backtrack if labels are consistent
            const Vector<Column> &labels = fetch_labels->get();
            Vector<Column> inter;
            label_set = &labels;
            std::set_intersection(label_intersection_.begin(), label_intersection_.end(),
                                  labels.begin(), labels.end(), std::back_inserter(inter));
            old_size = label_intersection_.size();
            cur_size = labels.size();
            new_size = inter.size();

            std::swap(label_intersection_, inter);

            if (label_intersection_.empty()) {
                cur_min_path_score_ = std::numeric_limits<score_t>::max();
                return;
            }
        }

        if (cur_size > new_size && this->prev_starts.count(trace[last_path_size_])) {
            // if the new coordinate or label set is a subset, then store
            // the difference to start backtracking from here later
            Vector<Column> diff;
            Vector<Tuple> diff_coords;
            auto prev_find = diff_label_sets_.find(trace[last_path_size_]);
            if (prev_find == diff_label_sets_.end()) {
                if (label_intersection_coords_.size()) {
                    utils::indexed_set_op<Tuple, CoordDiff>(
                        label_set->begin(), label_set->end(), label_coord->begin(),
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        std::back_inserter(diff), std::back_inserter(diff_coords)
                    );
                } else if (label_intersection_.size()) {
                    std::set_difference(label_set->begin(), label_set->end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                }
                if (diff.size()) {
                    diff_label_sets_.emplace(trace[last_path_size_],
                                             std::make_pair(std::move(diff),
                                                            std::move(diff_coords)));
                    this->prev_starts.erase(trace[last_path_size_]);
                }
            } else {
                auto &[prev_labels, prev_coords] = prev_find.value();
                if (label_intersection_coords_.size()) {
                    utils::indexed_set_op<Tuple, CoordDiff>(
                        prev_labels.begin(), prev_labels.end(), prev_coords.begin(),
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        std::back_inserter(diff), std::back_inserter(diff_coords)
                    );
                } else {
                    std::set_difference(prev_labels.begin(), prev_labels.end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                }
                std::swap(prev_labels, diff);
                std::swap(prev_coords, diff_coords);
                if (diff.size())
                    this->prev_starts.erase(trace[last_path_size_]);
            }
        }
    }

    cur_min_path_score_ = label_intersection_.size()
        ? extensions_.get_min_path_score(label_intersection_)
        : std::numeric_limits<score_t>::max();

    if ((!config_.allow_left_trim && table_i) || (
            config_.allow_left_trim && ((clipping && ops.data().back().first != Cigar::MATCH)
            || window.size() < this->config_.min_seed_length
            || alignment_score < std::max(min_path_score, cur_min_path_score_))))
        return;

    Alignment alignment = this->construct_alignment(
        ops, clipping, window, path, match, alignment_score, offset, extra_penalty
    );

    // store the label and coordinate information
    alignment.label_encoder = &labeled_graph_.get_annotator().get_label_encoder();

    alignment.label_columns = label_intersection_;

    if (labeled_graph_.get_coordinate_matrix()) {
        assert(label_intersection_.size() == label_intersection_coords_.size());
        alignment.label_coordinates = label_intersection_coords_;

        // if we were aligning backwards, then c represents the
        // end coordinate, so shift it
        if (coord_step == 1) {
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c -= path.size() - 1;
                }
            }
        }

        if (alignment.get_sequence().size() < this->graph_->get_k()) {
            size_t residual_offset = this->graph_->get_k() - alignment.get_sequence().size();
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c += residual_offset;
                }
            }
        }
    }

    if (!table_i && seed_labels_.size()) {
        Vector<Column> diff;
        if (seed_label_coordinates_.size()) {
            Vector<Tuple> diff_coords;
            utils::indexed_set_op<Tuple, CoordDiff>(
                seed_labels_.begin(), seed_labels_.end(), seed_label_coordinates_.begin(),
                alignment.label_columns.begin(), alignment.label_columns.end(),
                alignment.label_coordinates.begin(),
                std::back_inserter(diff), std::back_inserter(diff_coords)
            );
            std::swap(diff_coords, seed_label_coordinates_);
        } else {
            std::set_difference(seed_labels_.begin(), seed_labels_.end(),
                                alignment.label_columns.begin(),
                                alignment.label_columns.end(),
                                std::back_inserter(diff));
        }
        std::swap(diff, seed_labels_);
    }

    // instead of calling back the alignment, we add it to the local alignment buffer
    // after extension is done, the best ones are called back
    extensions_.add_alignment(std::move(alignment));
}

} // namespace align
} // namespace graph
} // namespace mtg
