#include "aligner_labeled.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/graph_extensions/node_lcs.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::binmat::BinaryMatrix::Row Row;
typedef annot::binmat::BinaryMatrix::Column Column;
typedef AnnotationBuffer::Columns Columns;
typedef DeBruijnGraph::node_index node_index;

// dummy index for an unfetched annotations
static constexpr size_t nannot = std::numeric_limits<size_t>::max();

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        column_sets_({ {} }) {
    if (multi_int_ && graph_.get_mode() != DeBruijnGraph::BASIC) {
        multi_int_ = nullptr;
        logger->warn("Coordinates not supported when aligning to CANONICAL "
                     "or PRIMARY mode graphs");
    }
}

void AnnotationBuffer::fetch_queued_annotations() {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    std::vector<node_index> queued_nodes;
    std::vector<Row> queued_rows;

    const DeBruijnGraph *base_graph = &graph_;

    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph_);
    if (canonical)
        base_graph = &canonical->get_graph();

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    for (const auto &path : queued_paths_) {
        std::vector<node_index> base_path;
        if (base_graph->get_mode() == DeBruijnGraph::CANONICAL) {
            // TODO: avoid this call of spell_path
            std::string query = spell_path(graph_, path);
            base_path = map_to_nodes(*base_graph, query);

        } else if (canonical) {
            base_path.reserve(path.size());
            for (node_index node : path) {
                base_path.emplace_back(canonical->get_base_node(node));
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
                if (node_to_cols_.emplace(path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                continue;
            }

            if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i]))) {
                // skip dummy nodes
                if (node_to_cols_.emplace(base_path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                        && base_path[i] != path[i]
                        && node_to_cols_.emplace(path[i], 0).second && has_coordinates()) {
                    label_coords_.emplace_back();
                }

                continue;
            }

            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            if (canonical || graph_.get_mode() == DeBruijnGraph::BASIC) {
                if (node_to_cols_.emplace(base_path[i], nannot).second) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                }

                continue;
            }

            assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);

            auto find_a = node_to_cols_.find(path[i]);
            auto find_b = node_to_cols_.find(base_path[i]);

            if (find_a == node_to_cols_.end() && find_b == node_to_cols_.end()) {
                node_to_cols_.emplace(path[i], nannot);
                queued_rows.push_back(row);
                queued_nodes.push_back(path[i]);

                if (path[i] != base_path[i]) {
                    node_to_cols_.emplace(base_path[i], nannot);
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                }
            } else if (find_a == node_to_cols_.end() && find_b != node_to_cols_.end()) {
                node_to_cols_.emplace(path[i], find_b->second);
                if (find_b->second == nannot) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(path[i]);
                }
            } else if (find_a != node_to_cols_.end() && find_b == node_to_cols_.end()) {
                node_to_cols_.emplace(base_path[i], find_a->second);
            } else {
                size_t label_i = std::min(find_a->second, find_b->second);
                if (label_i != nannot) {
                    find_a.value() = label_i;
                    find_b.value() = label_i;
                }
            }
        }
    }

    queued_paths_.clear();

    if (queued_nodes.empty())
        return;

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels) {
        assert(node_it != queued_nodes.end());
        assert(node_to_cols_.count(*node_it));
        assert(node_to_cols_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));

        size_t label_i = cache_column_set(std::move(labels));
        node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
        if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            assert(base_node == *node_it);
            node_to_cols_[*node_it] = label_i;
        } else if (canonical) {
            node_to_cols_[base_node] = label_i;
        } else {
            node_to_cols_[*node_it] = label_i;
            if (base_node != *node_it && node_to_cols_.emplace(base_node, label_i).second
                    && has_coordinates()) {
                label_coords_.emplace_back(label_coords_.back());
            }
        }
    };

    auto node_it = queued_nodes.begin();
    auto row_it = queued_rows.begin();
    if (has_coordinates()) {
        assert(multi_int_);
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows)) {
            Columns labels;
            labels.reserve(row_tuples.size());
            label_coords_.emplace_back();
            label_coords_.back().reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                label_coords_.back().push_back(std::move(coords));
            }
            push_node_labels(node_it++, row_it++, std::move(labels));
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
            push_node_labels(node_it++, row_it++, std::move(labels));
        }
    }

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols_) {
        assert(val != nannot);
    }
#endif
}

auto AnnotationBuffer::get_labels_and_coords(node_index node) const
        -> std::pair<const Columns*, const CoordinateSet*> {
    std::pair<const Columns*, const CoordinateSet*> ret_val { nullptr, nullptr };

    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph_))
        node = canonical->get_base_node(node);

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

template <class T1, class T2>
bool overlap_with_diff(const T1 &tuple1, const T2 &tuple2, size_t diff) {
    auto a_begin = tuple1.begin();
    const auto a_end = tuple1.end();
    auto b_begin = tuple2.begin();
    const auto b_end = tuple2.end();

    assert(std::is_sorted(a_begin, a_end));
    assert(std::is_sorted(b_begin, b_end));

    while (a_begin != a_end && b_begin != b_end) {
        if (*a_begin + diff == *b_begin)
            return true;

        if (*a_begin + diff < *b_begin) {
            ++a_begin;
        } else {
            ++b_begin;
        }
    }

    return false;
}

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

struct CoordIntersection {
    CoordIntersection(ssize_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + offset_ < *b_begin) {
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                *out = *a_begin;
                ++a_begin;
                ++b_begin;
                ++out;
            }
        }
    }

    ssize_t offset_;
};

LabeledExtender::LabeledExtender(const IDBGAligner &aligner, std::string_view query)
    : LabeledExtender(aligner.get_graph(), aligner.get_config(),
                      dynamic_cast<const LabeledAligner<>&>(aligner).annotation_buffer_,
                      query) {}

void LabeledExtender::flush() {
    annotation_buffer_.fetch_queued_annotations();
    for ( ; last_flushed_table_i_ < table.size(); ++last_flushed_table_i_) {
        auto &table_elem = table[last_flushed_table_i_];

        size_t parent_i = table_elem.parent_i;
        assert(parent_i < last_flushed_table_i_);

        auto clear = [&]() {
            DEBUG_LOG("Removed table element {}", last_flushed_table_i_);
            node_labels_[last_flushed_table_i_] = 0;
            node_labels_switched_[last_flushed_table_i_] = false;
            std::fill(table_elem.S.begin(), table_elem.S.end(), config_.ninf);
            std::fill(table_elem.E.begin(), table_elem.E.end(), config_.ninf);
            std::fill(table_elem.F.begin(), table_elem.F.end(), config_.ninf);
        };

        if (!node_labels_[parent_i]) {
            clear();
            continue;
        }

        if (node_labels_[parent_i] != node_labels_[last_flushed_table_i_]
                || node_labels_switched_[last_flushed_table_i_]) {
            continue;
        }

        node_index node = table_elem.node;
        const auto &parent_labels
            = annotation_buffer_.get_cached_column_set(node_labels_[parent_i]);

        auto cur_labels = annotation_buffer_.get_labels(node);
        assert(cur_labels);
        Columns intersect_labels;
        std::set_intersection(parent_labels.begin(), parent_labels.end(),
                              cur_labels->begin(), cur_labels->end(),
                              std::back_inserter(intersect_labels));
        if (intersect_labels.empty()) {
            clear();
        } else {
            node_labels_[last_flushed_table_i_]
                = annotation_buffer_.cache_column_set(std::move(intersect_labels));
        }
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (!DefaultColumnExtender::set_seed(seed))
        return false;

    assert(std::all_of(seed.get_nodes().begin(), seed.get_nodes().end(),
                       [&](node_index n) { return annotation_buffer_.get_labels(n); }));

    // the first node of the seed has already been flushed
    last_flushed_table_i_ = 1;

    remaining_labels_i_ = annotation_buffer_.cache_column_set(seed.label_columns);
    assert(remaining_labels_i_ != nannot);
    node_labels_.assign(1, remaining_labels_i_);
    node_labels_switched_.assign(1, false);
    base_coords_ = seed.label_coordinates;
    if (!base_coords_.size())
        return true;

    if (dynamic_cast<const RCDBG*>(graph_)) {
        for (auto &coords : base_coords_) {
            for (auto &c : coords) {
                c += seed.get_nodes().size() - 1 - seed.get_offset();
            }
        }
    } else if (seed.get_offset()) {
        for (auto &coords : base_coords_) {
            for (auto &c : coords) {
                c -= seed.get_offset();
            }
        }
    }

    return true;
}

template <class Outgoing, class LabelIntersectDiff>
double compute_label_change_scores(const DeBruijnGraph &graph,
                                   AnnotationBuffer &annotation_buffer,
                                   const DBGAlignerConfig &config,
                                   Alignment::node_index node,
                                   char node_last_char,
                                   const Outgoing &outgoing,
                                   LabelIntersectDiff &intersect_diff_labels) {
    assert(outgoing.size() == intersect_diff_labels.size());
    double sum_diff_probs = 0.0;

    size_t found_diffs = 0;
    for (const auto &[inter, diff, lclogprob] : intersect_diff_labels) {
        found_diffs += !diff.empty();
    }

    if (!found_diffs)
        return sum_diff_probs;

    auto default_init = [&]() {
        for (auto &[inter, diff, lclogprob] : intersect_diff_labels) {
            if (diff.size()) {
                sum_diff_probs += 1.0 - std::pow(2.0, config.label_change_score);
                lclogprob = std::floor(config.label_change_score - log2(outgoing.size()));
            }
        }
        sum_diff_probs = std::floor(log2(sum_diff_probs + 1.0) - log2(outgoing.size()));
    };

    // TODO: this cascade of graph unwrapping is ugly, find a cleaner way to do it
    const DeBruijnGraph *base_graph = &graph;
    if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(base_graph)) {
        base_graph = &rc_dbg->get_graph();
    } else if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_graph)) {
        base_graph = &canonical->get_graph();
    }

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);

    if (!dbg_succ) {
        default_init();
        return sum_diff_probs;
    }

    size_t expand_seed = graph.get_k() - 1
        - std::min(graph.get_k() - 1, config.label_change_search_width);

    std::shared_ptr<NodeLCS> lcs = dbg_succ->get_extension<NodeLCS>();

    if (!lcs) {
        default_init();
        return sum_diff_probs;
    }

    typedef boss::BOSS::TAlphabet TAlphabet;
    typedef boss::BOSS::edge_index edge_index;
    typedef Alignment::node_index node_index;
    typedef Alignment::Column Column;

    const boss::BOSS &boss = dbg_succ->get_boss();
    const auto *rev_graph = dynamic_cast<const RCDBG*>(&graph);
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(
        rev_graph ? &rev_graph->get_graph() : &graph
    );

    auto get_labeled_node = [&](node_index node) {
        if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
            return canonical->get_base_node(node);

        if (graph.get_mode() == DeBruijnGraph::CANONICAL) {
            std::string query = spell_path(graph, { node });
            return map_to_nodes(graph, query)[0];
        }

        return node;
    };

    assert(annotation_buffer.get_labels(node));
    node_index node_base = get_labeled_node(node);
    assert(node_base);

    // Case 1:  node == node_base, next == next_base
    //          node         next
    // e.g.,    GTCATGC C -> TCATGCC X
    // look for $$$ATGC C -> $$ATGCC X
    // Case 2:  node != node_base, next != next_base

    // Case 3:  node == node_base, next != next_base (TODO: handle later)
    // Case 4:  node != node_base, next == next_base (TODO: handle later)

    size_t num_nonbase_next_nodes = 0;
    for (const auto &[next, c, score] : outgoing) {
        num_nonbase_next_nodes += (next != get_labeled_node(next));
    }

    if (node != node_base && canonical && rev_graph && num_nonbase_next_nodes) {
        rev_graph = nullptr;
    } else {
        if (node != node_base) {
            default_init();
            return sum_diff_probs;
        }

        if (rev_graph && num_nonbase_next_nodes) {
            default_init();
            return sum_diff_probs;
        }
    }

    //          node         next
    // e.g.,    GTCATGC C <- XGTCATG C
    edge_index edge_base = dbg_succ->kmer_to_boss_index(
        rev_graph ? get_labeled_node(std::get<0>(outgoing[0])) : node_base
    );
    edge_index last_edge_base = boss.succ_last(edge_base);

    std::pair<edge_index, edge_index> node_range {
        boss.pred_last(last_edge_base - 1) + 1, last_edge_base
    };
    assert(node_range.first);
    assert(node_range.second);
    node_range = lcs->expand(node_range.first, node_range.second,
                             boss.get_k(), boss.get_k() - expand_seed);
    assert((*lcs)[node_range.first] < expand_seed);
    assert(node_range.second + 1 == lcs->data().size()
        || (*lcs)[node_range.second + 1] < expand_seed);

    auto next_range = node_range;
    TAlphabet cur_edge_label = rev_graph
        ? boss.get_node_last_value(dbg_succ->kmer_to_boss_index(node_base))
        : boss.encode(node_last_char);
    assert(cur_edge_label == boss.get_W(edge_base) % boss.alph_size);
    bool check = boss.tighten_range(&next_range.first, &next_range.second,
                                    cur_edge_label);
    std::ignore = check;
    assert(check);

    assert(boss.rank_W(node_range.second, cur_edge_label)
            - boss.rank_W(node_range.first - 1, cur_edge_label)
                == boss.rank_last(next_range.second) - boss.rank_last(next_range.first - 1));

    std::vector<TAlphabet> next_range_w;
    if (!rev_graph) {
        next_range_w.reserve(next_range.second - next_range.first + 1);
        for (edge_index i = next_range.first; i <= next_range.second; ++i) {
            next_range_w.emplace_back(boss.get_W(i) % boss.alph_size);
        }
    }

    std::vector<node_index> nodes_to_queue;
    for (edge_index i = boss.succ_W(node_range.first, cur_edge_label);
            i <= node_range.second;
            i = boss.succ_W(i + 1, cur_edge_label)) {
        if (node_index n = dbg_succ->boss_to_kmer_index(i))
            nodes_to_queue.push_back(n);
    }
    for (edge_index i = boss.succ_W(node_range.first, cur_edge_label + boss.alph_size);
            i <= node_range.second;
            i = boss.succ_W(i + 1, cur_edge_label + boss.alph_size)) {
        if (node_index n = dbg_succ->boss_to_kmer_index(i))
            nodes_to_queue.push_back(n);
    }

    for (edge_index i = next_range.first; i <= next_range.second; ++i) {
        if (node_index n = dbg_succ->boss_to_kmer_index(i))
            nodes_to_queue.push_back(n);
    }

    annotation_buffer.queue_path(std::move(nodes_to_queue));
    annotation_buffer.fetch_queued_annotations();

    size_t total = 0;
    Vector<size_t> out_map(boss.alph_size, intersect_diff_labels.size());
    for (size_t i = 0; i < outgoing.size(); ++i) {
        const auto &[next, c, score] = outgoing[i];
        out_map[boss.encode(c)] = i;
    }

    Vector<size_t> counts(boss.alph_size);
    Vector<tsl::hopscotch_set<Column>> exclude_sets;
    exclude_sets.reserve(intersect_diff_labels.size());
    for (const auto &[inter, diff, _] : intersect_diff_labels) {
        exclude_sets.emplace_back(diff.begin(), diff.end());
        exclude_sets.back().insert(inter.begin(), inter.end());
    }

    edge_index node_wp = boss.succ_W(node_range.first, cur_edge_label);
    edge_index node_wm = boss.succ_W(node_range.first, cur_edge_label + boss.alph_size);
    edge_index i = boss.succ_last(next_range.first);
    edge_index i_start = boss.pred_last(i - 1) + 1;
    bool label_preserving_traversal_found = false;
    auto step = [&]() {
        if (node_wp < node_wm) {
            node_wp = boss.succ_W(node_wp + 1, cur_edge_label);
        } else {
            node_wm = boss.succ_W(node_wm + 1, cur_edge_label + boss.alph_size);
        }

        if (node_wp <= node_wm) {
            i_start = i + 1;
            i = boss.succ_last(i_start);
        }
    };
    for ( ; i <= next_range.second; step()) {
        assert(node_wp != node_wm);
        edge_index node_w = std::min(node_wp, node_wm);
        assert(boss.fwd(node_w, cur_edge_label) == i);
        node_index n = dbg_succ->boss_to_kmer_index(node_w);
        if (!n)
            continue;

        ++total;
        sdsl::bit_vector match_found(boss.alph_size, false);
        size_t matches = 0;
        const Alignment::Columns *n_labels = annotation_buffer.get_labels(n);
        assert(n_labels);

        // TODO: find a more efficient way to get this
        TAlphabet common_next_c_enc = rev_graph
            ? boss.encode(complement(boss.decode(
                boss.get_minus_k_value(node_w, boss.get_k() - 1).first)
              ))
            : 0;

        for (edge_index next_edge = i_start; next_edge <= i; ++next_edge) {
            if (matches == match_found.size() - 1)
                break;

            node_index m = dbg_succ->boss_to_kmer_index(next_edge);
            if (!m)
                continue;

            TAlphabet next_c_enc = rev_graph
                ? common_next_c_enc
                : next_range_w[next_edge - next_range.first];
            size_t j = out_map[next_c_enc];
            if (j == outgoing.size() || match_found[next_c_enc])
                continue;

            const auto &[_, diff, lclogprob] = intersect_diff_labels[j];
            if (diff.empty())
                continue;

            const auto &exclude = exclude_sets[j];
            const Alignment::Columns *m_labels = annotation_buffer.get_labels(m);
            assert(m_labels);
            Alignment::Columns inter;
            std::set_intersection(n_labels->begin(), n_labels->end(),
                                  m_labels->begin(), m_labels->end(),
                                  std::back_inserter(inter));
            for (Column col : inter) {
                if (!exclude.count(col)) {
                    // found an independent label
                    match_found[next_c_enc] = true;
                    ++matches;
                }
            }
        }

        if (matches)
            label_preserving_traversal_found = true;

        for (size_t j = 0; j < match_found.size(); ++j) {
            counts[j] += match_found[j];
        }
    }

    if (!label_preserving_traversal_found)
        return sum_diff_probs;

    for (size_t i = 0; i < outgoing.size(); ++i) {
        const auto &[next, c, score] = outgoing[i];
        auto &[inter, diff, lclogprob] = intersect_diff_labels[i];
        if (size_t count = counts[boss.encode(c)]) {
            assert(diff.size());
            lclogprob = log2(count) - log2(total);
            sum_diff_probs += 1.0 - std::pow(2.0, lclogprob);
            lclogprob = std::floor(lclogprob - log2(outgoing.size()));
        }
    }

    return std::floor(log2(sum_diff_probs + 1.0) - log2(outgoing.size()));
}

void LabeledExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    assert(table.size() == node_labels_.size());
    assert(table.size() == node_labels_switched_.size());

    // if we are in the seed and want to force the seed to be fixed, automatically
    // take the next node in the seed
    size_t next_offset = table[table_i].offset + 1;
    bool is_seed_fixed = force_fixed_seed || fixed_seed();
    bool in_seed = next_offset - seed_->get_offset() < seed_->get_sequence().size()
                    && (next_offset < graph_->get_k() || is_seed_fixed);

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
        [&](node_index next, char c, score_t score) {
            outgoing.emplace_back(next, c, score);
            if (!in_seed)
                annotation_buffer_.queue_path({ next });
        },
        table_i, force_fixed_seed
    );

    if (outgoing.empty())
        return;

    if (outgoing.size() == 1) {
        const auto &[next, c, score] = outgoing[0];
        if (next_offset <= graph_->get_k()
                || next_offset - graph_->get_k() - 1 >= seed_->label_column_diffs.size()
                || seed_->label_column_diffs[next_offset - graph_->get_k() - 1].empty()) {
            // Assume that annotations are preserved in unitigs. Violations of this
            // assumption are corrected after the next flush
            node_labels_.emplace_back(node_labels_[table_i]);
            node_labels_switched_.emplace_back(config_.label_change_union);
        } else {
            const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);
            const auto &diff_labels = seed_->label_column_diffs[next_offset - graph_->get_k() - 1];
            Columns decoded;
            std::set_symmetric_difference(columns.begin(), columns.end(),
                                          diff_labels.begin(), diff_labels.end(),
                                          std::back_inserter(decoded));
            node_labels_.emplace_back(annotation_buffer_.cache_column_set(std::move(decoded)));
            node_labels_switched_.emplace_back(true);
            last_flushed_table_i_ = std::max(last_flushed_table_i_, node_labels_.size());
        }

        callback(next, c, score);
        return;
    }

    // flush the AnnotationBuffer and correct for annotation errors introduced above
    flush();

    if (!node_labels_[table_i])
        return;

    assert(annotation_buffer_.get_labels(node));

    // use the label set of the current node in the alignment tree as the basis
    const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);

    // no coordinates are present in the annotation
    if (!annotation_buffer_.get_labels_and_coords(node).second) {
        std::vector<std::tuple<Columns, Columns, double>> intersect_diff_labels(
            outgoing.size()
        );

        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (size_t i = 0; i < outgoing.size(); ++i) {
            const auto &[next, c, score] = outgoing[i];
            auto &[intersect_labels, diff_labels, lclogprob] = intersect_diff_labels[i];
            lclogprob = config_.ninf;
            auto next_labels = annotation_buffer_.get_labels(next);
            assert(next_labels);

            set_intersection_difference(next_labels->begin(), next_labels->end(),
                                        columns.begin(), columns.end(),
                                        std::back_inserter(intersect_labels),
                                        std::back_inserter(diff_labels));

            if (!config_.label_change_union && intersect_labels.size())
                diff_labels.clear();
        }

        double sum_diff_probs = compute_label_change_scores(
            *graph_, annotation_buffer_, config_, node, table[table_i].c, outgoing,
            intersect_diff_labels
        );

        for (size_t i = 0; i < outgoing.size(); ++i) {
            const auto &[next, c, score] = outgoing[i];
            score_t abs_match_score = std::abs(config_.score_matrix[c][c]);
            auto &[inter, diff, lclogprob] = intersect_diff_labels[i];
            assert(diff.size() || lclogprob == config_.ninf);

            if (lclogprob > config_.ninf && diff.size())
                lclogprob *= abs_match_score;

            DEBUG_LOG("Position {} edge {}: label preserve: {}\tlabel change: {}",
                      next_offset - seed_->get_offset(), i, sum_diff_probs * abs_match_score,
                      lclogprob);

            if (config_.label_change_union) {
                if (lclogprob > config_.ninf && diff.size()) {
                    Columns column_union;
                    std::set_union(columns.begin(), columns.end(), diff.begin(), diff.end(),
                                   std::back_inserter(column_union));
                    node_labels_.emplace_back(annotation_buffer_.cache_column_set(std::move(column_union)));
                    node_labels_switched_.emplace_back(true);
                    callback(next, c, score + lclogprob);
                } else {
                    node_labels_.emplace_back(node_labels_[table_i]);
                    node_labels_switched_.emplace_back(true);
                    callback(next, c, score + sum_diff_probs * abs_match_score);
                }
            } else {
                if (lclogprob > config_.ninf && diff.size()) {
                    node_labels_.emplace_back(annotation_buffer_.cache_column_set(std::move(diff)));
                    node_labels_switched_.emplace_back(true);
                    callback(next, c, score + lclogprob);
                } else if (inter.size()) {
                    node_labels_.emplace_back(annotation_buffer_.cache_column_set(std::move(inter)));
                    node_labels_switched_.emplace_back(false);
                    callback(next, c, score + sum_diff_probs * abs_match_score);
                }
            }
        }

        return;
    }

    // check label and coordinate consistency
    // use the seed as the basis for labels and coordinates
    assert(seed_->label_coordinates.size());

    // compute the coordinate distance from base_coords
    size_t dist = next_offset - graph_->get_k() + 1;

    for (const auto &[next, c, score] : outgoing) {
        const Columns *base_labels = &seed_->label_columns;
        const CoordinateSet *base_coords = &base_coords_;
        auto [next_labels, next_coords]
            = annotation_buffer_.get_labels_and_coords(next);

        assert(next_coords);

        // if we are traversing backwards, then negate the coordinate delta
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(base_labels, next_labels);
            std::swap(base_coords, next_coords);
        }

        // check if at least one label has consistent coordinates
        Columns intersect_labels;

        try {
            auto col_it = columns.begin();
            auto col_end = columns.end();

            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                next_labels->begin(), next_labels->end(), next_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    // also, intersect with the label set of node
                    while (col_it != col_end && c > *col_it) {
                        ++col_it;
                    }

                    if (col_it == col_end)
                        throw std::exception();

                    if (c < *col_it)
                        return;

                    // then check coordinate consistency with the seed
                    if (overlap_with_diff(coords, other_coords, dist))
                        intersect_labels.push_back(c);
                }
            );
        } catch (const std::exception&) {}

        if (intersect_labels.size()) {
            // found a consistent set of coordinates, so assign labels for this next node
            node_labels_.push_back(annotation_buffer_.cache_column_set(
                                                std::move(intersect_labels)));
            node_labels_switched_.emplace_back(false);
            callback(next, c, score);
        }
    }
}

bool LabeledExtender::skip_backtrack_start(size_t i) {
    assert(remaining_labels_i_ != nannot);
    assert(node_labels_[i] != nannot);

    // if this alignment tree node has been visited previously, ignore it
    assert(remaining_labels_i_);
    if (!prev_starts.emplace(i).second)
        return true;

    if (std::any_of(node_labels_switched_.begin(), node_labels_switched_.end(),
                    [](const auto &a) { return a; })) {
        return false;
    }

    // check if this starting point involves seed labels which have not been considered yet
    const auto &end_labels = annotation_buffer_.get_cached_column_set(node_labels_[i]);
    const auto &left_labels = annotation_buffer_.get_cached_column_set(remaining_labels_i_);

    label_intersection_ = Columns{};
    label_diff_ = Columns{};
    set_intersection_difference(left_labels.begin(), left_labels.end(),
                                end_labels.begin(), end_labels.end(),
                                std::back_inserter(label_intersection_),
                                std::back_inserter(label_diff_));
    label_diff_.push_back(nannot);

    return label_intersection_.empty();
}

void LabeledExtender::call_alignments(score_t end_score,
                                      const std::vector<node_index> &path,
                                      const std::vector<size_t> &trace,
                                      const std::vector<score_t> &score_trace,
                                      const Cigar &ops,
                                      size_t clipping,
                                      size_t offset,
                                      std::string_view window,
                                      const std::string &match,
                                      score_t extra_score,
                                      const std::function<void(Alignment&&)> &callback) {
    assert(label_intersection_.size()
            || std::any_of(node_labels_switched_.begin(), node_labels_switched_.end(),
                           [](const auto &a) { return a; }));

    Alignment alignment = construct_alignment(ops, clipping, window, path, match,
                                              end_score, offset, score_trace, extra_score);
    alignment.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();

    auto [base_labels, base_coords]
        = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().front());
    assert(base_labels);
    assert(base_labels->size());

    if (!clipping)
        base_labels = &seed_->label_columns;

    auto call_alignment = [&]() {
        assert(alignment.label_columns.size());
        if (label_diff_.size() && label_diff_.back() == nannot) {
            label_diff_.pop_back();
            remaining_labels_i_ = annotation_buffer_.cache_column_set(std::move(label_diff_));
            assert(remaining_labels_i_ != nannot);
            label_diff_ = Columns{};
        }

        callback(std::move(alignment));
    };

    if (!base_coords) {
        if (std::any_of(node_labels_switched_.begin(), node_labels_switched_.end(),
                        [](const auto &a) { return a; })) {
            auto it = trace.rend() - alignment.get_nodes().size();
            assert(*it);
            alignment.label_columns
                = annotation_buffer_.get_cached_column_set(node_labels_[*it]);
            assert(alignment.label_columns.size());

            alignment.label_column_diffs.reserve(alignment.get_nodes().size() - 1);
            const auto *prev = &annotation_buffer_.get_cached_column_set(node_labels_[*it]);
            for (++it; it != trace.rend(); ++it) {
                const auto *cur = &annotation_buffer_.get_cached_column_set(node_labels_[*it]);
                alignment.label_column_diffs.emplace_back();
                std::set_symmetric_difference(
                    prev->begin(), prev->end(), cur->begin(), cur->end(),
                    std::back_inserter(alignment.label_column_diffs.back())
                );
                std::swap(prev, cur);
            }

            if (alignment.extra_scores.empty())
                alignment.extra_scores.resize(alignment.label_column_diffs.size(), 0);

            label_diff_.assign(1, nannot);
        } else {
            alignment.label_columns = std::move(label_intersection_);
            label_intersection_ = Columns{};
        }

        call_alignment();
        return;
    }

    ssize_t dist = alignment.get_nodes().size() - 1;
    if (!clipping) {
        base_coords = &seed_->label_coordinates;
        dist -= seed_->get_offset();
        if (dynamic_cast<const RCDBG*>(graph_))
            dist = alignment.get_sequence().size() - seed_->get_sequence().size();
    }

    auto label_it = label_intersection_.begin();
    auto label_end_it = label_intersection_.end();

    if (alignment.get_nodes().size() == 1) {
        auto it = base_labels->begin();
        auto end = base_labels->end();
        auto c_it = base_coords->begin();
        while (label_it != label_end_it && it != end) {
            if (*label_it < *it) {
                ++label_it;
            } else if (*label_it > *it) {
                ++it;
                ++c_it;
            } else {
                alignment.label_columns.emplace_back(*it);
                alignment.label_coordinates.emplace_back(*c_it);
                ++it;
                ++c_it;
                ++label_it;
            }
        }
    } else {
        auto [cur_labels, cur_coords]
            = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().back());
        assert(cur_labels);
        assert(cur_labels->size());
        assert(cur_coords);
        assert(cur_coords->size());
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(cur_labels, base_labels);
            std::swap(cur_coords, base_coords);
        }
        CoordIntersection intersect_coords(dist);
        try {
            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                cur_labels->begin(), cur_labels->end(), cur_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    while (label_it != label_end_it && c > *label_it) {
                        ++label_it;
                    }

                    if (label_it == label_end_it)
                        throw std::exception();

                    if (c < *label_it)
                        return;

                    Alignment::Tuple overlap;
                    intersect_coords(coords.begin(), coords.end(),
                                     other_coords.begin(), other_coords.end(),
                                     std::back_inserter(overlap));
                    if (overlap.size()) {
                        alignment.label_columns.emplace_back(c);
                        alignment.label_coordinates.emplace_back(std::move(overlap));
                    }
                }
            );
        } catch (const std::exception&) {}
    }

    if (alignment.label_coordinates.empty())
        return;

    if (dynamic_cast<const RCDBG*>(graph_) && alignment.get_offset()) {
        for (auto &coords : alignment.label_coordinates) {
            for (auto &c : coords) {
                c += alignment.get_offset();
            }
        }
    }

    call_alignment();
}

template <class Seeder, class Extender, class AlignmentCompare>
LabeledAligner<Seeder, Extender, AlignmentCompare>
::LabeledAligner(const DeBruijnGraph &graph,
                 const DBGAlignerConfig &config,
                 const Annotator &annotator)
      : DBGAligner<Seeder, Extender, AlignmentCompare>(graph, config),
        annotation_buffer_(graph, annotator),
        annotator_(annotator) {
    // do not use a global xdrop cutoff since we need separate cutoffs for each label
    if (annotation_buffer_.has_coordinates())
        this->config_.global_xdrop = false;
}

template <class Seeder, class Extender, class AlignmentCompare>
LabeledAligner<Seeder, Extender, AlignmentCompare>::~LabeledAligner() {
    logger->trace("Buffered {}/{} nodes and {} label combinations",
                  annotation_buffer_.num_nodes_buffered(),
                  this->graph_.num_nodes(),
                  annotation_buffer_.num_column_sets());
}

template <class Seeder, class Extender, class AlignmentCompare>
auto LabeledAligner<Seeder, Extender, AlignmentCompare>
::build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                const std::vector<AlignmentResults> &wrapped_seqs) const -> BatchSeeders {
    BatchSeeders seeders
        = DBGAligner<Seeder, Extender, AlignmentCompare>::build_seeders(seq_batch, wrapped_seqs);

    // now we're going to filter the seeds
    logger->trace("Filtering seeds by label");
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds;
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds_rc;

    size_t num_seeds = 0;
    size_t num_seeds_rc = 0;

    for (auto &[seeder, nodes, seeder_rc, nodes_rc] : seeders) {
        counted_seeds.emplace_back(seeder->get_seeds(), seeder->get_num_matches());
        num_seeds += counted_seeds.back().first.size();

        auto add_seeds = [&](const auto &seeds) {
            for (const Alignment &seed : seeds) {
                annotation_buffer_.queue_path(std::vector<node_index>(seed.get_nodes()));
            }
        };

        add_seeds(counted_seeds.back().first);

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            counted_seeds_rc.emplace_back(seeder_rc->get_seeds(),
                                          seeder_rc->get_num_matches());
            num_seeds_rc += counted_seeds_rc.back().first.size();
            add_seeds(counted_seeds_rc.back().first);
        }
#endif
    }

    logger->trace("Prefetching labels");
    annotation_buffer_.fetch_queued_annotations();

    size_t num_seeds_left = 0;
    size_t num_seeds_rc_left = 0;

    for (size_t i = 0; i < counted_seeds.size(); ++i) {
        auto &[seeder, nodes, seeder_rc, nodes_rc] = seeders[i];
        auto &[seeds, num_matching] = counted_seeds[i];
        if (seeds.size()) {
            num_matching = filter_seeds(seeds);
            num_seeds_left += seeds.size();
        }

        seeder = make_shared<ManualSeeder>(std::move(seeds), num_matching);

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            assert(seeder_rc);
            auto &[seeds, num_matching] = counted_seeds_rc[i];
            if (seeds.size()) {
                num_matching = filter_seeds(seeds);
                num_seeds_rc_left += seeds.size();
            }

            seeder_rc = make_shared<ManualSeeder>(std::move(seeds), num_matching);
        }
#endif
    }

    logger->trace("Kept {}/{} seeds", num_seeds_left + num_seeds_rc_left,
                                      num_seeds + num_seeds_rc);

    return seeders;
}

inline size_t get_num_matches(const std::vector<Alignment> &seeds) {
    size_t num_matching = 0;
    size_t last_end = 0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty())
            continue;

        size_t begin = seeds[i].get_clipping();
        size_t end = begin + seeds[i].get_query_view().size();
        if (end > last_end) {
            num_matching += end - begin;
            if (begin < last_end)
                num_matching -= last_end - begin;
        }

        if (size_t offset = seeds[i].get_offset()) {
            size_t clipping = seeds[i].get_clipping();
            for (++i; i < seeds.size()
                        && seeds[i].get_offset() == offset
                        && seeds[i].get_clipping() == clipping; ++i) {}
            --i;
        }

        last_end = end;
    }
    return num_matching;
}

template <class Seeder, class Extender, class AlignmentCompare>
size_t LabeledAligner<Seeder, Extender, AlignmentCompare>
::filter_seeds(std::vector<Alignment> &seeds) const {
    VectorMap<Column, uint64_t> label_counter;
    for (const Alignment &seed : seeds) {
        for (node_index node : seed.get_nodes()) {
            assert(annotation_buffer_.get_labels(node));
            for (uint64_t label : *annotation_buffer_.get_labels(node)) {
                ++label_counter[label];
            }
        }
    }

    if (label_counter.empty())
        return get_num_matches(seeds);

    std::vector<std::pair<Column, uint64_t>> label_counts
        = const_cast<std::vector<std::pair<Column, uint64_t>>&&>(
            label_counter.values_container()
        );

    std::sort(label_counts.begin(), label_counts.end(), utils::GreaterSecond());

    double cutoff = static_cast<double>(label_counts[0].second)
                        * this->config_.rel_score_cutoff;
    auto it = std::find_if(label_counts.begin() + 1, label_counts.end(),
                           [cutoff](const auto &a) { return a.second < cutoff; });

    label_counts.erase(it, label_counts.end());

    Columns labels;
    labels.reserve(label_counts.size());
    for (const auto &[label, count] : label_counts) {
        labels.push_back(label);
    }
    std::sort(labels.begin(), labels.end());

    for (Alignment &seed : seeds) {
        const std::vector<node_index> &nodes = seed.get_nodes();
        auto [fetch_labels, fetch_coords]
            = annotation_buffer_.get_labels_and_coords(nodes[0]);
        if (!fetch_labels)
            continue;

        if (fetch_coords) {
            auto a_begin = fetch_labels->begin();
            auto a_end = fetch_labels->end();
            auto a_c_begin = fetch_coords->begin();
            auto b_begin = labels.begin();
            auto b_end = labels.end();
            while (a_begin != a_end && b_begin != b_end) {
                if (*a_begin < *b_begin) {
                    ++a_begin;
                    ++a_c_begin;
                } else if (*a_begin > *b_begin) {
                    ++b_begin;
                } else {
                    seed.label_columns.emplace_back(*a_begin);
                    seed.label_coordinates.emplace_back(*a_c_begin);
                    ++a_begin;
                    ++a_c_begin;
                    ++b_begin;
                }
            }
        } else {
            std::set_intersection(fetch_labels->begin(), fetch_labels->end(),
                                  labels.begin(), labels.end(),
                                  std::back_inserter(seed.label_columns));
        }

        seed.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();

        for (size_t i = 1; i < nodes.size() && seed.label_columns.size(); ++i) {
            auto [next_fetch_labels, next_fetch_coords]
                = annotation_buffer_.get_labels_and_coords(nodes[i]);
            if (next_fetch_coords) {
                Columns label_inter;
                Alignment::CoordinateSet coord_inter;
                CoordIntersection intersect_coords(i);
                utils::match_indexed_values(
                    seed.label_columns.begin(), seed.label_columns.end(),
                    seed.label_coordinates.begin(),
                    next_fetch_labels->begin(), next_fetch_labels->end(),
                    next_fetch_coords->begin(),
                    [&](auto col, const auto &coords, const auto &other_coords) {
                        Alignment::Tuple overlap;
                        intersect_coords(coords.begin(), coords.end(),
                                         other_coords.begin(), other_coords.end(),
                                         std::back_inserter(overlap));
                        if (overlap.size()) {
                            label_inter.push_back(col);
                            coord_inter.push_back(std::move(overlap));
                        }
                    }
                );

                std::swap(seed.label_columns, label_inter);
                std::swap(seed.label_coordinates, coord_inter);
            } else if (next_fetch_labels) {
                Columns temp;
                std::set_intersection(next_fetch_labels->begin(), next_fetch_labels->end(),
                                      seed.label_columns.begin(), seed.label_columns.end(),
                                      std::back_inserter(temp));
                std::swap(temp, seed.label_columns);
            } else {
                seed.label_columns.clear();
                seed.label_coordinates.clear();
            }
        }

        annotation_buffer_.cache_column_set(seed.label_columns);

        if (seed.get_offset() && seed.label_coordinates.size()) {
            for (auto &tuple : seed.label_coordinates) {
                for (auto &coord : tuple) {
                    coord += seed.get_offset();
                }
            }
        }
    }

    auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [&](const auto &a) {
        return a.label_columns.empty();
    });

    seeds.erase(seed_it, seeds.end());

    return get_num_matches(seeds);
}

template class LabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
