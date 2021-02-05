#include "aligner_dp_table.hpp"

#include <tsl/hopscotch_set.h>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType>
bool DPTable<NodeType>::add_seed(const Alignment<NodeType> &seed,
                                 const DBGAlignerConfig &config,
                                 size_t size,
                                 size_t start_pos,
                                 size_t query_offset) {
    query_offset_ = query_offset;
    char start_char = *(seed.get_query_end() - 1);
    score_t last_char_score = config.get_row(start_char)[seed.get_sequence().back()];

    iterator column_it = dp_table_.find(seed.back());
    if (column_it == dp_table_.end()) {
        column_it = emplace(seed.back(), Column(size, config.min_cell_score,
                                                start_char, start_pos)).first;
    } else {
        expand_to_cover(column_it, 0, size);
    }

    auto &table_init = column_it.value();

    bool update = false;

    if (table_init.best_score() < seed.get_score()) {
        auto last_op = seed.get_cigar().back().first;
        table_init.scores[start_pos] = seed.get_score();
        table_init.ops[start_pos] = last_op;
        table_init.prev_nodes[start_pos] = 0;
        table_init.gap_prev_nodes[start_pos] = 0;

        if (last_op != Cigar::DELETION && last_op != Cigar::MATCH
                && last_op != Cigar::MISMATCH) {
            throw std::runtime_error("Seeds must end in DELETION, MATCH, or MISMATCH");
        }

        // This vector stores the best scores for partial alignments ending in a
        // delete operation. If the last operation in the seed is MATCH or MISMATCH,
        // then we replace it with a score corresponding to INSERTION followed by DELETION
        table_init.gap_scores[start_pos] = std::max(
            last_op == Cigar::DELETION
                ? table_init.scores[start_pos]
                : table_init.scores[start_pos] - last_char_score
                                               + config.gap_opening_penalty
                                               + config.gap_opening_penalty,
            config.min_cell_score
        );

        table_init.gap_count[start_pos] = 1;
        update = true;
    }

    return update;
}

template <typename NodeType>
void DPTable<NodeType>
::extract_alignments(const DeBruijnGraph &graph,
                     const DBGAlignerConfig &config,
                     std::string_view query_view,
                     std::function<void(Alignment<NodeType>&&, NodeType)> callback,
                     score_t min_path_score,
                     const Alignment<NodeType> &seed,
                     NodeType *node) {
    NodeType start_node;
    if (config.num_alternative_paths == 1 && node) {
        // avoid sorting column iterators if we're only interested in the top path
        start_node = DeBruijnGraph::npos;
        auto column_it = dp_table_.find(*node);
        assert(column_it != dp_table_.end());
        Alignment<NodeType> alignment(*this,
                                      config,
                                      query_view,
                                      column_it,
                                      column_it->second.best_pos,
                                      graph.get_k() - 1,
                                      &start_node,
                                      seed);

        if (alignment.empty() && !alignment.get_query().data())
            return;

        assert(alignment.is_valid(graph, &config));

        callback(std::move(alignment), start_node);

        return;
    }

    // store visited nodes in paths to avoid returning subalignments
    tsl::hopscotch_set<key_type> visited_nodes;

    std::vector<const_iterator> starts;
    starts.reserve(size());
    for (const_iterator it = dp_table_.begin(); it != dp_table_.end(); ++it) {
        if (it->second.best_score() > min_path_score
                && it->second.best_op() == Cigar::MATCH)
            starts.emplace_back(it);
    }

    if (starts.empty())
        return;

    std::sort(starts.begin(), starts.end(),
              [](const auto &a, const auto &b) {
                  return a->second.best_score() > b->second.best_score();
              });

    size_t num_paths = 0;
    for (const_iterator column_it : starts) {
        if (num_paths >= config.num_alternative_paths)
            return;

        // ignore if the current point is a subalignment of one already visited
        if (visited_nodes.find(column_it->first) != visited_nodes.end())
            continue;

        start_node = DeBruijnGraph::npos;
        Alignment<NodeType> next(*this,
                                 config,
                                 query_view,
                                 column_it,
                                 column_it->second.best_pos,
                                 graph.get_k() - 1,
                                 &start_node,
                                 seed);

        if (next.empty() && !next.get_query().data())
            continue;

        assert(next.is_valid(graph, &config));
        visited_nodes.insert(next.begin(), next.end());

        callback(std::move(next), start_node);

        ++num_paths;
    }
}

template <typename NodeType>
void DPTable<NodeType>::Column::expand_to_cover(size_t begin, size_t end) {
    assert(best_pos >= start_index);
    assert(best_pos - start_index < scores.size());
    assert(last_priority_pos >= start_index);
    assert(last_priority_pos - start_index < scores.size());

    if (begin >= start_index) {
        // the current range already covers [begin, end)
        if (end <= start_index + scores.size() - 8)
            return;

        // extend the range to the right to reach end
        scores.resize(end + 8 - start_index, min_score_);
        gap_scores.resize(end + 8 - start_index, min_score_);
        ops.resize(end + 8 - start_index, Cigar::CLIPPED);
        prev_nodes.resize(end + 8 - start_index, 0);
        gap_prev_nodes.resize(end + 8 - start_index, 0);
        gap_count.resize(end + 8 - start_index, 0);
    } else if (end <= start_index + scores.size() - 8) {
        // extend the range to the left to reach begin
        size_t shift = start_index - begin;
        start_index = begin;
        scores.insert(scores.begin(), shift, min_score_);
        gap_scores.insert(gap_scores.begin(), shift, min_score_);
        ops.insert(ops.begin(), shift, Cigar::CLIPPED);
        prev_nodes.insert(prev_nodes.begin(), shift, 0);
        gap_prev_nodes.insert(gap_prev_nodes.begin(), shift, 0);
        gap_count.insert(gap_count.begin(), shift, 0);
    } else {
        // extend the range in both directions
        size_t shift = start_index - begin;
        start_index = begin;
        end += 8;

        size_t new_size = end - begin;

        scores.reserve(new_size);
        gap_scores.reserve(new_size);
        ops.reserve(new_size);
        prev_nodes.reserve(new_size);
        gap_prev_nodes.reserve(new_size);
        gap_count.reserve(new_size);

        scores.insert(scores.begin(), shift, min_score_);
        gap_scores.insert(gap_scores.begin(), shift, min_score_);
        ops.insert(ops.begin(), shift, Cigar::CLIPPED);
        prev_nodes.insert(prev_nodes.begin(), shift, 0);
        gap_prev_nodes.insert(gap_prev_nodes.begin(), shift, 0);
        gap_count.insert(gap_count.begin(), shift, 0);

        scores.resize(new_size, min_score_);
        gap_scores.resize(new_size, min_score_);
        ops.resize(new_size, Cigar::CLIPPED);
        prev_nodes.resize(new_size, 0);
        gap_prev_nodes.resize(new_size, 0);
        gap_count.resize(new_size, 0);
    }

    assert(best_pos >= start_index);
    assert(best_pos - start_index < scores.size());
    assert(last_priority_pos >= start_index);
    assert(last_priority_pos - start_index < scores.size());
}

template class DPTable<>;

} // namespace align
} // namespace graph
} // namespace mtg
