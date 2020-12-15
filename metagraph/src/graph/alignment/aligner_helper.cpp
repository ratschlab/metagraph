#include "aligner_helper.hpp"

#include <tsl/hopscotch_set.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;


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
        table_init.gap_scores[start_pos] = std::max(
            last_op == Cigar::INSERTION
                ? table_init.scores[start_pos]
                : table_init.scores[start_pos] - last_char_score + config.gap_opening_penalty + config.gap_opening_penalty,
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
                     const std::string_view query_view,
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
Alignment<NodeType>::Alignment(const std::string_view query,
                               std::vector<NodeType>&& nodes,
                               std::string&& sequence,
                               score_t score,
                               size_t clipping,
                               bool orientation,
                               size_t offset)
      : query_begin_(query.data()),
        query_end_(query.data() + query.size()),
        nodes_(std::move(nodes)),
        sequence_(std::move(sequence)),
        score_(score),
        orientation_(orientation),
        offset_(offset) {
    assert(query_begin_);
    assert(query_end_);

    size_t query_size = query_end_ - query_begin_;
    size_t min_length = std::min(query_size, sequence_.size());

    cigar_ = std::inner_product(
        query_begin_,
        query_begin_ + min_length,
        sequence_.c_str(),
        Cigar(Cigar::CLIPPED, clipping),
        [&](Cigar &cigar, bool equal) -> Cigar& {
            cigar.append(equal ? Cigar::MATCH : Cigar::MISMATCH);
            return cigar;
        },
        std::equal_to<char>()
    );

    assert(!(query_size - min_length) || (sequence_.size() - min_length));
    cigar_.append(Cigar::DELETION, query_size - min_length);
    cigar_.append(Cigar::INSERTION, sequence_.size() - min_length);
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const DPTable<NodeType> &dp_table,
                               const DBGAlignerConfig &config,
                               const std::string_view query_view,
                               typename DPTable<NodeType>::const_iterator column,
                               size_t start_pos,
                               size_t offset,
                               NodeType *start_node,
                               const Alignment &seed)
      : query_begin_(NULL),
        query_end_(NULL),
        orientation_(seed.get_orientation()),
        offset_(offset) {
    assert(start_node);

    auto i = start_pos;
    size_t shift = column->second.start_index;
    assert(i >= shift);
    assert(i - shift < column->second.scores.size());
    score_ = column->second.scores.at(i - shift);
    Cigar::Operator op = column->second.ops.at(i - shift);
    NodeType prev_node;
    switch (column->second.prev_nodes.at(i - shift)) {
        case 0: { prev_node = SequenceGraph::npos; } break;
        case 0xFF: { prev_node = column->first; } break;
        default: {
            prev_node = column->second.select_prev_node(column->second.prev_nodes.at(i - shift));
        }
    }

    NodeType prev_gap_node;
    switch (column->second.gap_prev_nodes.at(i - shift)) {
        case 0: { prev_gap_node = SequenceGraph::npos; } break;
        case 0xFF: { prev_gap_node = column->first; } break;
        default: {
            prev_gap_node = column->second.select_prev_node(column->second.gap_prev_nodes.at(i - shift));
        }
    }

    if (op == Cigar::INSERTION)
        prev_node = prev_gap_node;

    uint32_t gap_count = op == Cigar::INSERTION ? column->second.gap_count.at(i - shift) - 1 : 0;

    if (!i && prev_node == SequenceGraph::npos)
        return;

    // use config to recompute CIGAR score in DEBUG mode
    score_t score_track = score_;
    Cigar::Operator last_op = Cigar::CLIPPED;

    score_t gap_diff = config.gap_opening_penalty - config.gap_extension_penalty;

    // TODO: If there is a cyclic part of the graph in which the optimal
    //       alignment involves a deletion, then a score which was previously
    //       from a deletion may be replaced with a match. This will cause
    //       subsequent deletion scores to be wrong since they're no longer
    //       extensions. The only way to fix this is to store a separate vector
    //       to keep partial alignments ending in deletions.
    //       Until this is fixed, the score checking asserts have been commented out.

    std::vector<typename DPTable<NodeType>::const_iterator> out_columns;
    while (prev_node != SequenceGraph::npos) {
        auto prev_column = dp_table.find(prev_node);
        assert(prev_column != dp_table.end());
        assert(i || op == Cigar::INSERTION);

        switch (op) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                --i;
                out_columns.emplace_back(column);

                if (last_op == Cigar::INSERTION)
                    score_track -= gap_diff;

                // assert(column->second.scores.at(i + 1) >= score_track);
                score_track -= config.get_row(column->second.last_char)[query_view[i]];
                // assert(prev_column->second.scores.at(i) >= score_track);

            } break;
            case Cigar::INSERTION: {
                out_columns.emplace_back(column);

                score_track -= config.gap_extension_penalty;

            } break;
            case Cigar::DELETION: {
                assert(column == prev_column);
                --i;

                assert(i >= shift);
                assert(i - shift < column->second.scores.size());

                // assert(column->second.prev_nodes.at(i + 1) == 0xFF);
                // assert(column->second.scores.at(i + 1) >= score_track);
                assert(column->second.ops.at(i - shift) != Cigar::INSERTION);
                score_track -= column->second.ops.at(i - shift) == Cigar::DELETION
                    ? config.gap_extension_penalty
                    : config.gap_opening_penalty;
                // assert(column->second.scores.at(i) >= score_track);

            } break;
            case Cigar::CLIPPED: { assert(false); }
        }

        cigar_.append(op);

        last_op = op;

        column = prev_column;
        shift = prev_column->second.start_index;
        assert(i >= shift);
        assert(i - shift < column->second.scores.size());
        if (gap_count) {
            --gap_count;
        } else {
            op = column->second.ops.at(i - shift);

            if (op == Cigar::INSERTION)
                gap_count = column->second.gap_count.at(i - shift) - 1;
        }
        switch (column->second.prev_nodes.at(i - shift)) {
            case 0: { prev_node = SequenceGraph::npos; } break;
            case 0xFF: { prev_node = column->first; } break;
            default: {
                prev_node = column->second.select_prev_node(column->second.prev_nodes.at(i - shift));
            }
        }

        switch (column->second.gap_prev_nodes.at(i - shift)) {
            case 0: { prev_gap_node = SequenceGraph::npos; } break;
            case 0xFF: { prev_gap_node = column->first; } break;
            default: {
                prev_gap_node = column->second.select_prev_node(column->second.gap_prev_nodes.at(i - shift));
            }
        }

        if (op == Cigar::INSERTION)
            prev_node = prev_gap_node;
    }

    const auto &score_col = column->second.scores;

    if (last_op == Cigar::INSERTION)
        score_track -= gap_diff;

    score_t correction = score_col.at(i - shift) - score_track;

    // assert(correction >= 0);

    if (correction > 0)
        logger->trace("Fixing outdated score: {} -> {}", score_, score_ + correction);

    score_ -= score_col.at(i - shift) - correction;

    *start_node = column->first;

    if (i > std::numeric_limits<Cigar::LengthType>::max())
        throw std::runtime_error("Error: clipping length can't be stored in CIGAR");

    cigar_.append(Cigar::CLIPPED, i);
    assert(cigar_.size());

    query_begin_ = query_view.data() + i;
    query_end_ = query_view.data() + start_pos;

    std::reverse(cigar_.begin(), cigar_.end());


    nodes_.resize(out_columns.size());
    std::transform(out_columns.rbegin(),
                   out_columns.rend(),
                   nodes_.begin(),
                   [](const auto &iter) { return iter->first; });

    sequence_.assign(out_columns.size(), 'N');
    std::transform(out_columns.rbegin(),
                   out_columns.rend(),
                   sequence_.begin(),
                   [](const auto &iter) { return iter->second.last_char; });

    if (correction < 0) {
        logger->warn("Incorrect score found: {} -> {}\nQuery: {}\nTarget: {}\nCIGAR: {}",
                     score_, score_ + correction,
                     seed.get_sequence() + std::string(query_view),
                     seed.get_sequence() + sequence_,
                     seed.get_cigar().to_string() + cigar_.to_string());
    }
}

template <typename NodeType>
void Alignment<NodeType>::append(Alignment&& other) {
    assert(query_end_ == other.query_begin_);
    assert(orientation_ == other.orientation_);
    assert(cigar_.empty() || cigar_.back().first != Cigar::CLIPPED);

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;

    cigar_.append(std::move(other.cigar_));
    query_end_ = other.query_end_;
}


template <typename NodeType>
AlignmentPrefix<NodeType>& AlignmentPrefix<NodeType>::operator++() {
    assert(cigar_it_ < cigar_rend_);

    switch (*cigar_it_) {
        case Cigar::Operator::MATCH:
        case Cigar::Operator::MISMATCH: {
            assert(end_it_ > begin_it_);
            assert(ref_end_it_ > ref_begin_it_);
            assert(node_it_ > alignment_->begin());
            score_ -= config_->get_row(*(end_it_ - 1))[*(ref_end_it_ - 1)];
            --end_it_;
            --ref_end_it_;

            // if we're at the first node of the alignment, then traverse
            // backwards to reconstruct a path
            if (node_it_ - 1 == alignment_->begin()) {
                auto last_node = prefix_node_ ? prefix_node_ : alignment_->front();
                assert(last_node);
                prefix_node_ = DeBruijnGraph::npos;
                graph_->adjacent_incoming_nodes(last_node, [&](auto prev) {
                    // TODO: what if there are multiple?
                    if (!prefix_node_)
                        prefix_node_ = prev;
                });
                ++offset_;
            } else {
                --node_it_;
            }
        } break;
        case Cigar::Operator::DELETION: {
            assert(end_it_ > begin_it_);
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

            --end_it_;
        } break;
        case Cigar::Operator::INSERTION: {
            assert((--(cigar_it_.base())).get_it() >= alignment_->get_cigar().begin());
            assert(ref_end_it_ > ref_begin_it_);
            assert(node_it_ > alignment_->begin());
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

            --ref_end_it_;

            // if we're at the first node of the alignment, then traverse
            // backwards to reconstruct a path
            if (node_it_ - 1 == alignment_->begin()) {
                auto last_node = prefix_node_ ? prefix_node_ : alignment_->front();
                assert(last_node);
                prefix_node_ = DeBruijnGraph::npos;
                graph_->adjacent_incoming_nodes(last_node, [&](auto prev) {
                    // TODO: what if there are multiple?
                    if (!prefix_node_)
                        prefix_node_ = prev;
                });
                ++offset_;

            } else {
                --node_it_;
            }
        } break;
        case Cigar::Operator::CLIPPED: { assert(false); }
    }

    ++cigar_it_;
    ++trim_;

    assert(cigar_rbegin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_rend_);

    return *this;
}

template <typename NodeType>
AlignmentPrefix<NodeType>& AlignmentPrefix<NodeType>::operator--() {
    throw std::runtime_error("not implemented");
    return *this;
}

template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator++() {
    assert(cigar_it_ < cigar_end_);
    switch (*cigar_it_) {
        case Cigar::Operator::MATCH:
        case Cigar::Operator::MISMATCH: {
            assert(begin_it_ < end_it_);
            assert(ref_begin_it_ < ref_end_it_);
            score_ -= config_->get_row(*begin_it_)[*ref_begin_it_];
            ++begin_it_;
            ++ref_begin_it_;
        } break;
        case Cigar::Operator::DELETION: {
            assert(begin_it_ < end_it_);
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            ++begin_it_;
        } break;
        case Cigar::Operator::INSERTION: {
            assert(ref_begin_it_ < ref_end_it_);
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            ++ref_begin_it_;
        } break;
        case Cigar::Operator::CLIPPED: { assert(false); }
    }

    ++cigar_it_;
    ++trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    return *this;
}

template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator--() {
    assert(cigar_it_.get_it() >= alignment_->get_cigar().begin());

    --cigar_it_;

    switch (*cigar_it_) {
        case Cigar::Operator::MATCH:
        case Cigar::Operator::MISMATCH: {
            assert(begin_it_ >= alignment_->get_query().data());
            assert(ref_begin_it_ >= alignment_->get_sequence().data());
            --begin_it_;
            --ref_begin_it_;
            score_ += config_->get_row(*begin_it_)[*ref_begin_it_];
        } break;
        case Cigar::Operator::DELETION: {
            assert(begin_it_ >= alignment_->get_query().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            --begin_it_;
        } break;
        case Cigar::Operator::INSERTION: {
            assert(cigar_it_.get_it() >= alignment_->get_cigar().begin());
            assert(ref_begin_it_ >= alignment_->get_sequence().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            --ref_begin_it_;
        } break;
        case Cigar::Operator::CLIPPED: { assert(false); }
    }

    --trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    return *this;
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const AlignmentPrefix<NodeType> &alignment_prefix) {
    const auto &data = alignment_prefix.data();
    assert(data.get_cigar().size());

    offset_ = data.get_offset() + alignment_prefix.get_offset();
    if (offset_ == alignment_prefix.get_graph().get_k()) {
        *this = Alignment();
        return;
    }

    query_begin_ = data.get_query().data();
    auto prefix_node = alignment_prefix.get_prefix_node();
    if (prefix_node) {
        nodes_.assign(&prefix_node, &prefix_node + 1);
    } else {
        nodes_.assign(data.get_nodes().begin(), alignment_prefix.get_node_end_it());
    }

    cigar_ = data.get_cigar();
    orientation_ = data.get_orientation();
    sequence_ = alignment_prefix.get_sequence();

    score_ = alignment_prefix.get_score();
    query_end_ = alignment_prefix.get_query().data() + alignment_prefix.get_query().size();

    if (cigar_.back().first == Cigar::CLIPPED)
        cigar_.pop_back();

    for (size_t i = 0; i < alignment_prefix.get_trim(); ++i) {
        --cigar_.back().second;

        if (!cigar_.back().second)
            cigar_.pop_back();
    }
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const AlignmentSuffix<NodeType> &alignment_suffix) {
    const auto &data = alignment_suffix.data();
    assert(data.get_cigar().size());

    query_end_ = data.get_query().data() + data.get_query().size();
    orientation_ = data.get_orientation();
    offset_ = data.get_offset();

    score_ = alignment_suffix.get_score();
    query_begin_ = alignment_suffix.get_query().data();
    sequence_ = alignment_suffix.get_sequence();

    auto node_it = data.begin();
    const auto node_end_it = data.end();

    CigarOpIterator begin(data.get_cigar(), data.get_clipping());
    assert(begin <= alignment_suffix.get_op_it());
    CigarOpIterator end(data.get_cigar(),
                        data.get_cigar().end() - static_cast<bool>(data.get_end_clipping()));
    assert(begin <= end);
    assert(alignment_suffix.get_op_it() <= end);

    auto it = alignment_suffix.get_op_it();
    assert(begin <= it);
    assert(it <= end);

    while (begin < it) {
        if (*begin != Cigar::DELETION) {
            if (node_it + 1 != node_end_it) {
                ++node_it;
                if (offset_)
                    --offset_;
            } else {
                ++offset_;
            }
        }
        ++begin;
    }

    if (offset_ == alignment_suffix.get_k()) {
        *this = Alignment();
        return;
    }

    while (it != end) {
        assert(it < end);
        cigar_.append(*it);
        ++it;
    }

    nodes_.assign(node_it, node_end_it);

    extend_query_begin(data.get_query().data());
}

template <typename NodeType>
std::pair<Alignment<NodeType>, Alignment<NodeType>> Alignment<NodeType>
::get_best_overlap(const Alignment &first, const Alignment &second,
                   const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config) {
    assert(first.is_valid(graph, &config));
    assert(second.is_valid(graph, &config));
    assert(first.sequence_.size() - first.nodes_.size() + 1 + first.offset_
        == second.sequence_.size() - second.nodes_.size() + 1 + second.offset_);

    // no overlap
    if (second.get_query().data() >= first.get_query_end())
        return std::make_pair(first, second);

    // check if first (second) is a prefix of second (first)
    if (first.get_query().data() == second.get_query().data()) {
        auto [first_it, second_it] = std::mismatch(
            CigarOpIterator(first.cigar_),
            CigarOpIterator(first.cigar_, first.cigar_.end()),
            CigarOpIterator(second.cigar_)
        );
        if (first_it.get_it() == first.cigar_.end()) {
            // first cigar is a prefix of the second
            return std::make_pair(Alignment(), second);
        } else if (second_it.get_it() == second.cigar_.end()) {
            // second cigar is a prefix of the first
            return std::make_pair(first, Alignment());
        }
    }

    // check if first (second) is a suffix of second (first)
    if (first.get_query_end() == second.get_query_end()) {
        auto first_rbegin = std::make_reverse_iterator(CigarOpIterator(first.cigar_, first.cigar_.end()));
        auto first_rend = std::make_reverse_iterator(CigarOpIterator(first.cigar_, first.get_clipping()));
        auto second_rbegin = std::make_reverse_iterator(CigarOpIterator(second.cigar_, second.cigar_.end()));
        auto second_rend = std::make_reverse_iterator(CigarOpIterator(second.cigar_, second.get_clipping()));

        auto [first_rit, second_rit] = std::mismatch(first_rbegin, first_rend, second_rbegin);
        if (first_rit == first_rend) {
            // first cigar is a suffix of the second
            return std::make_pair(Alignment(), second);
        } else if (second_rit == second_rend) {
            // second cigar is a suffix of the first
            return std::make_pair(first, Alignment());
        }
    }

    // Now, find the maximal scoring combination of a prefix of the first alignment
    // and a suffix of the second alignment such that they do not overlap.
    // When trimming, INSERTION only consumes the reference, so these can be discarded

    AlignmentPrefix<NodeType> first_prefix(first, config, graph);
    AlignmentPrefix<NodeType> best_first_prefix = first_prefix;
    AlignmentSuffix<NodeType> second_suffix(second, config, graph.get_k());
    AlignmentSuffix<NodeType> best_second_suffix = second_suffix;

    size_t overlap = first_prefix.get_query().data() + first_prefix.get_query().size()
        - second_suffix.get_query().data();

    score_t best_score = config.min_cell_score;
    score_t cur_score;
    size_t best_overlap = overlap + 1;

    // initialize: trim overlap from the left of second
    for (size_t i = 0; i < overlap; ++i) {
        assert(!second_suffix.eof());
        while (second_suffix.get_front_op() == Cigar::Operator::INSERTION) {
            assert(!second_suffix.eof());
            ++second_suffix;
            assert(Alignment(second_suffix).is_valid(graph, &config));
        }
        assert(!second_suffix.eof());
        ++second_suffix;
        assert(Alignment(second_suffix).is_valid(graph, &config));
        while (!second_suffix.eof() && second_suffix.get_front_op() == Cigar::Operator::INSERTION) {
            assert(!second_suffix.eof());
            ++second_suffix;
            assert(Alignment(second_suffix).is_valid(graph, &config));
        }
        assert(!second_suffix.eof() || i + 1 == overlap);
    }

    if (second_suffix.get_query().size() >= config.min_seed_length) {
        cur_score = first_prefix.get_score() + second_suffix.get_score();
        if (cur_score > best_score) {
            best_first_prefix = first_prefix;
            best_second_suffix = second_suffix;
            best_score = cur_score;
            best_overlap = 0;
        }
    }

    // iteratively shift the splice point
    auto node_it = first.begin();
    for (size_t i = 1; i <= overlap; ++i) {
        assert(first_prefix.get_query().data() + first_prefix.get_query().size()
            == second_suffix.get_query().data());

        while (!first_prefix.eof()
                && first_prefix.get_back_op() == Cigar::Operator::INSERTION
                && first_prefix.get_node_end_it() > node_it) {
            ++first_prefix;
            assert(Alignment(first_prefix).is_valid(graph, &config));
            cur_score = first_prefix.get_score() + second_suffix.get_score();
            if (cur_score > best_score || (cur_score >= best_score && best_overlap == i)) {
                best_first_prefix = first_prefix;
                best_second_suffix = second_suffix;
                best_score = cur_score;
                best_overlap = i;
            }
        }

        if (first_prefix.eof() || first_prefix.get_node_end_it() == node_it)
            break;

        ++first_prefix;
        assert(Alignment(first_prefix).is_valid(graph, &config));

        while (!first_prefix.eof()
                && first_prefix.get_back_op() == Cigar::Operator::INSERTION
                && first_prefix.get_node_end_it() > node_it) {
            ++first_prefix;
            assert(Alignment(first_prefix).is_valid(graph, &config));
        }

        if (!second_suffix.eof()) {
            while (!second_suffix.reof()
                    && second_suffix.get_front_op() == Cigar::Operator::INSERTION) {
                --second_suffix;
                assert(Alignment(second_suffix).is_valid(graph, &config));
            }
        }

        assert(!second_suffix.reof());

        --second_suffix;
        assert(Alignment(second_suffix).is_valid(graph, &config));

        while (!second_suffix.reof()
                && second_suffix.get_front_op() == Cigar::Operator::INSERTION) {
            --second_suffix;
            assert(Alignment(second_suffix).is_valid(graph, &config));
        }

        if (first_prefix.get_query().size() < config.min_seed_length)
            break;

        if (second_suffix.get_query().size() >= config.min_seed_length) {
            cur_score = first_prefix.get_score() + second_suffix.get_score();
            if (cur_score > best_score) {
                best_first_prefix = first_prefix;
                best_second_suffix = second_suffix;
                best_score = cur_score;
                best_overlap = i;
            }
        }
    }

    assert(first_prefix.get_query().data() + first_prefix.get_query().size()
            == second_suffix.get_query().data());

    if (best_overlap > overlap)
        return std::make_pair(Alignment(), Alignment());

    if (best_second_suffix.get_query().data() == second.get_query_end())
        return std::make_pair(first, Alignment());

    if (best_first_prefix.get_query().data() + best_first_prefix.get_query().size()
            <= first.get_query().data()) {
        return std::make_pair(Alignment(), second);
    }

    Alignment<NodeType> new_first(best_first_prefix);
    Alignment<NodeType> new_second(best_second_suffix);

    new_second.extend_query_begin(second.get_query().data() - second.get_clipping());

    assert(new_first.get_query_end() == new_second.get_query().data());

    return std::make_pair(std::move(new_first), std::move(new_second));
}

template <typename NodeType>
void Alignment<NodeType>::trim_offset() {
    if (!offset_ || empty() || cigar_.empty())
        return;

    auto it = cigar_.begin();
    if (it->first == Cigar::CLIPPED)
        ++it;

    if (it == cigar_.end())
        return;

    auto jt = nodes_.begin();
    size_t counter = 0;
    while (offset_ && it != cigar_.end() && jt != nodes_.end()) {
        if (counter == it->second
                || it->first == Cigar::CLIPPED
                || it->first == Cigar::DELETION) {
            ++it;
            counter = 0;
            continue;
        }

        size_t jump = std::min(std::min(offset_, size_t(it->second)),
                               static_cast<size_t>(nodes_.end() - jt));
        offset_ -= jump;
        counter += jump;
        jt += jump;
    }

    if (jt == nodes_.end()) {
        --jt;
        ++offset_;
    }

    nodes_.erase(nodes_.begin(), jt);
}

template <typename NodeType>
void Alignment<NodeType>::reverse_complement(const DeBruijnGraph &graph,
                                             const std::string_view query_rev_comp) {
    assert(graph.is_canonical_mode());

    if (empty())
        return;

    assert(query_end_ + get_end_clipping()
        == query_begin_ - get_clipping() + query_rev_comp.size());
    assert(is_valid(graph));

    trim_offset();
    assert(is_valid(graph));

    if (!offset_) {
        reverse_complement_seq_path(graph, sequence_, nodes_);
    } else {
        assert(nodes_.size() == 1);
        // extract target sequence prefix
        std::string rev_seq = graph.get_node_sequence(nodes_.front()).substr(0, offset_)
            + sequence_;

        // if the alignment starts from a source k-mer, then traverse forwards
        // until a non-dummy k-mer is hit and check if its reverse complement exists
        const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
        if (dbg_succ && rev_seq[0] == boss::BOSS::kSentinel) {
            const auto &boss = dbg_succ->get_boss();
            assert(rev_seq.length() == graph.get_k());
            auto edge = dbg_succ->kmer_to_boss_index(nodes_.back());
            auto edge_label = boss.get_W(edge) % boss.alph_size;

            // TODO: for efficiency, the last outgoing edge is taken at each step.
            // This is valid for canonical graphs since the reverse complement of
            // the found non-dummy k-mer is guaranteed to exist, but this node may
            // not exist in a non-canonical graph.
            for (size_t i = 0; i < offset_; ++i) {
                edge = boss.fwd(edge, edge_label);
                edge_label = boss.get_W(edge) % boss.alph_size;
                if (edge_label == boss::BOSS::kSentinelCode) {
                    // reverse complement not found
                    *this = Alignment();
                    return;
                }

                nodes_.push_back(dbg_succ->boss_to_kmer_index(edge));
                rev_seq.push_back(boss.decode(edge_label));
            }
            nodes_.assign(nodes_.begin() + offset_, nodes_.end());
            rev_seq.assign(rev_seq.begin() + offset_, rev_seq.end());

            assert(nodes_ == map_sequence_to_nodes(graph, rev_seq));
            std::vector<NodeType> rev_nodes = nodes_;
            reverse_complement_seq_path(graph, rev_seq, rev_nodes);

            assert(std::find(rev_nodes.begin(), rev_nodes.end(),
                             DeBruijnGraph::npos) == rev_nodes.end());

            assert(rev_seq.size() > offset_);
            sequence_.assign(rev_seq.begin() + offset_, rev_seq.end());
            nodes_.assign(rev_nodes.begin(), rev_nodes.end());

        } else {
            assert(nodes_.size() == 1);
            assert(nodes_ == map_sequence_to_nodes(graph, rev_seq));
            std::vector<NodeType> rev_nodes = nodes_;
            reverse_complement_seq_path(graph, rev_seq, rev_nodes);

            assert(std::find(rev_nodes.begin(), rev_nodes.end(),
                             DeBruijnGraph::npos) == rev_nodes.end());

            // trim off ending from reverse complement (corresponding to the added prefix)
            size_t trim_left = offset_;
            while (trim_left && rev_nodes.size() > 1) {
                rev_nodes.pop_back();
                rev_seq.pop_back();
                --trim_left;
            }

            for (size_t i = 0; i < trim_left; ++i) {
                size_t indegree = 0;
                graph.adjacent_incoming_nodes(rev_nodes[0], [&](NodeType prev) {
                    ++indegree;

                    // TODO: there are multiple possible reverse complements, which
                    // do we pick? Currently we pick the first one
                    if (indegree == 1)
                        rev_nodes[0] = prev;
                });

                if (!indegree) {
                    *this = Alignment();
                    return;
                }

                rev_seq.pop_back();
            }

            nodes_ = rev_nodes;
            sequence_ = rev_seq;
            offset_ = trim_left;
            assert(!trim_left
                    || graph.get_node_sequence(rev_nodes[0]).substr(trim_left) == sequence_);
        }
    }

    std::reverse(cigar_.begin(), cigar_.end());

    orientation_ = !orientation_;

    query_begin_ = query_rev_comp.data() + get_clipping();
    query_end_ = query_rev_comp.data() + (query_rev_comp.size() - get_end_clipping());

    assert(query_end_ >= query_begin_);
    assert(is_valid(graph));
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
template <typename NodeType>
Json::Value Alignment<NodeType>::path_json(size_t node_size,
                                           const std::string &label) const {
    assert(nodes_.size());

    Json::Value path;

    auto cigar_it = cigar_.begin();
    if (cigar_.size() && cigar_it->first == Cigar::CLIPPED) {
        cigar_it++;
    }

    size_t cigar_offset = 0;
    assert(cigar_it != cigar_.end());

    int64_t rank = 1;
    auto query_start = query_begin_;

    size_t cur_pos = rank == 1 ? offset_ : 0;

    Json::Value mapping;
    Json::Value position;
    position["node_id"] = nodes_.front();

    if (cur_pos)
        position["offset"] = Json::Value::UInt64(cur_pos);

    // set to true if the node is the reverse complement of the query
    //position["is_reverse"] = false;

    mapping["position"] = position;

    // handle alignment to the first node
    while (cur_pos < node_size && cigar_it != cigar_.end()) {
        assert(cigar_it->second > cigar_offset);
        size_t next_pos = std::min(node_size,
                                   cur_pos + (cigar_it->second - cigar_offset));
        size_t next_size = next_pos - cur_pos;
        assert(cigar_offset + next_size <= cigar_it->second);

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::MISMATCH: {
                assert(query_start + next_size <= query_end_);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::DELETION: {
                assert(query_start + next_size <= query_end_);
                // this assumes that DELETIONS can't happen right after insertions
                //edit["from_length"] = 0;
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;

                // the target is not consumed, so reset the position
                next_pos = cur_pos;
            } break;
            case Cigar::INSERTION: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                assert(query_start + next_size <= query_end_);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                query_start += next_size;
            } break;
            case Cigar::CLIPPED: {
                ++cigar_it;
                cigar_offset = 0;
                assert(cigar_it == cigar_.end());
                continue;
            }
        }

        cigar_offset += next_size;
        cur_pos = next_pos;

        mapping["edit"].append(edit);

        if (cigar_offset == cigar_it->second) {
            ++cigar_it;
            cigar_offset = 0;
        }
    }

    mapping["rank"] = rank++;
    path["mapping"].append(mapping);

    // handle the rest of the alignment
    for (auto node_it = nodes_.begin() + 1; node_it != nodes_.end(); ++node_it) {
        assert(cigar_it != cigar_.end());
        assert(cigar_it->second > cigar_offset);

        Json::Value mapping;
        Json::Value position;
        position["node_id"] = *node_it;
        position["offset"] = Json::Value::UInt64(node_size - 1);
        // set to true if the node is the reverse complement of the query
        //position["is_reverse"] = false;
        mapping["position"] = position;

        if (cigar_it->first == Cigar::DELETION) {
            Json::Value edit;
            size_t length = cigar_it->second - cigar_offset;
            assert(query_start + length < query_end_);
            // TODO: this assumes that DELETIONS can't happen right after insertions
            //edit["from_length"] = 0;
            edit["to_length"] = Json::Value::UInt64(length);
            edit["sequence"] = std::string(query_start, length);
            query_start += length;
            ++cigar_it;
            cigar_offset = 0;
            mapping["edit"].append(edit);
            assert(cigar_it != cigar_.end());
        }

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::MISMATCH: {
                assert(query_start < query_end_);
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                edit["sequence"] = std::string(query_start, 1);
                query_start++;
            } break;
            case Cigar::INSERTION: {
                edit["from_length"] = 1;
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                query_start++;
            } break;
            case Cigar::DELETION:
            case Cigar::CLIPPED: assert(false); break;
        }

        if (++cigar_offset == cigar_it->second) {
            cigar_offset = 0;
            ++cigar_it;
        }

        mapping["edit"].append(edit);
        mapping["rank"] = rank++;
        path["mapping"].append(mapping);
    }

    assert(query_start == query_end_);
    assert(cigar_it == cigar_.end()
            || (cigar_it + 1 == cigar_.end() && cigar_it->first == Cigar::CLIPPED));

    path["length"] = Json::Value::UInt64(nodes_.size());
    //path["is_circular"]; // bool

    if (label.size())
        path["name"] = label;

    return path;
}

template <typename NodeType>
Json::Value Alignment<NodeType>::to_json(const std::string &query,
                                         const DeBruijnGraph &graph,
                                         bool is_secondary,
                                         const std::string &read_name,
                                         const std::string &label) const {
    assert(is_valid(graph));

    // encode alignment
    Json::Value alignment;

    alignment["name"] = read_name;
    alignment["sequence"] = query;

    if (sequence_.size())
        alignment["annotation"]["ref_sequence"] = sequence_;

    if (query_end_ == query_begin_)
        return alignment;

    auto query_start = query.c_str();
    assert(query_start + cigar_.get_clipping() == query_begin_);
    assert(query_end_ >= query_start);

    alignment["annotation"]["cigar"] = cigar_.to_string();

    // encode path
    if (nodes_.size())
        alignment["path"] = path_json(graph.get_k(), label);

    alignment["score"] = static_cast<int32_t>(score_);

    if (query_begin_ > query_start)
        alignment["query_position"] = static_cast<int32_t>(query_begin_ - query_start);

    if (is_secondary)
        alignment["is_secondary"] = is_secondary;

    alignment["identity"] = query_end_ != query_begin_
        ? static_cast<double>(get_num_matches()) / (query_end_ - query_begin_)
        : 0;

    alignment["read_mapped"] = (query_end_ != query_begin_);

    if (orientation_)
        alignment["read_on_reverse_strand"] = orientation_;


    if (cigar_.get_clipping())
        alignment["soft_clipped"] = static_cast<bool>(cigar_.get_clipping());

    // Unused flags (for now)
    //alignment["quality"]; // bytes
    //alignment["mapping_quality"]; // int32
    //alignment["sample_name"]; // string
    //alignment["read_group"]; // string
    //alignment["fragment_prev"]; // Alignment
    //alignment["fragment_next"]; // Alignment
    //alignment["fragment"]; // Path
    //alignment["refpos"]; // Position
    //alignment["read_paired"]; // bool
    //alignment["mate_unmapped"]; // bool
    //alignment["mate_on_reverse_strand"]; // bool
    //alignment["discordant_insert_size"]; // bool
    //alignment["uniqueness"]; // double
    //alignment["correct"]; // double
    //alignment["secondary_score"]; // repeated int32
    //alignment["fragment_score"]; // double
    //alignment["mate_mapped_to_disjoint_subgraph"]; // bool
    //alignment["fragment_length_distribution"]; // string
    //alignment["haplotype_scored"]; // bool
    //alignment["haplotype_logprob"]; // double
    //alignment["time_used"]; // double
    //alignment["to_correct"]; // Position
    //alignment["correctly_mapped"]; // bool

    return alignment;
}

template <typename NodeType>
std::shared_ptr<const std::string> Alignment<NodeType>
::load_from_json(const Json::Value &alignment,
                 const DeBruijnGraph &graph) {
    cigar_.clear();
    nodes_.clear();
    sequence_.clear();

    auto query_sequence = std::make_shared<const std::string>(
        alignment["sequence"].asString()
    );
    auto query_start = query_sequence->c_str();

    query_begin_ = query_start + alignment["query_position"].asInt();
    orientation_ = alignment["read_on_reverse_strand"].asBool();
    score_ = alignment["score"].asInt();

    const auto& mapping = alignment["path"]["mapping"];
    assert(mapping.size() == alignment["path"]["length"].asUInt64());

    Json::ArrayIndex i = 0;
    offset_ = mapping[i]["position"]["offset"].asUInt64();

    size_t path_steps = 0;

    query_end_ = query_begin_;

    cigar_.append(Cigar::CLIPPED, query_begin_ - query_start);
    for (; i < mapping.size(); ++i) {
        nodes_.emplace_back(mapping[i]["position"]["node_id"].asUInt64());
        if (nodes_.size() == 1) {
            sequence_ = graph.get_node_sequence(nodes_.back()).substr(offset_);
        } else {
            graph.call_outgoing_kmers(
                *(nodes_.rbegin() + 1),
                [&](auto node, char c) {
                    if (node == nodes_.back())
                        sequence_ += c;
                }
            );
        }
        const auto& edits = mapping[i]["edit"];

        for (Json::ArrayIndex j = 0; j < edits.size(); ++j) {
            assert(edits[j]["from_length"].asUInt64()
                       || edits[j]["to_length"].asUInt64());

            if (edits[j]["from_length"] == edits[j]["to_length"]) {
                if (edits[j]["sequence"].asString().size()) {
                    cigar_.append(Cigar::MISMATCH, edits[j]["from_length"].asUInt64());
                } else {
                    cigar_.append(Cigar::MATCH, edits[j]["from_length"].asUInt64());
                }

                path_steps += edits[j]["from_length"].asUInt64();
                query_end_ += edits[j]["to_length"].asUInt64();
            } else if (edits[j]["from_length"].asUInt64()) {
                cigar_.append(Cigar::INSERTION, edits[j]["from_length"].asUInt64());
                path_steps += edits[j]["from_length"].asUInt64();
            } else {
                assert(edits[j]["to_length"].asUInt64());
                cigar_.append(Cigar::DELETION, edits[j]["to_length"].asUInt64());
                query_end_ += edits[j]["to_length"].asUInt64();
            }
        }
    }

    if (query_end_ < query_sequence->c_str() + query_sequence->length()) {
        cigar_.append(Cigar::CLIPPED,
                      query_sequence->c_str() + query_sequence->length() - query_end_);
    }

    if (alignment["annotation"]["cigar"]) {
        assert(alignment["annotation"]["cigar"].asString() == cigar_.to_string());

        if (cigar_ != Cigar(alignment["annotation"]["cigar"].asString()))
            throw std::runtime_error("ERROR: CIGAR and mapping mismatch");
    }

    sequence_ = sequence_.substr(0, path_steps);
    assert(!alignment["annotation"]["ref_sequence"]
        || alignment["annotation"]["ref_sequence"].asString() == sequence_);

    if (!is_valid(graph))
        throw std::runtime_error("ERROR: JSON reconstructs invalid alignment");

    return query_sequence;
}

template <typename NodeType>
bool spell_path(const DeBruijnGraph &graph,
                const std::vector<NodeType> &path,
                std::string &seq,
                size_t offset = 0) {
    assert(offset < graph.get_k());

    if (path.empty())
        return "";

    if (std::find(path.begin(), path.end(), DeBruijnGraph::npos) != path.end()) {
        std::cerr << "ERROR: path has invalid nodes\n";

        for (NodeType node : path) {
            std::cerr << node << " ";
        }

        std::cerr << std::endl;

        return false;
    }

    seq.clear();
    seq.reserve(path.size() + graph.get_k() - 1 - offset);

    seq += graph.get_node_sequence(path.front()).substr(offset);

    for (size_t i = 1; i < path.size(); ++i) {
        char next = '\0';
        graph.call_outgoing_kmers(path[i - 1], [&](auto next_node, char c) {
            if (next_node == path[i])
                next = c;
        });

        if (!next) {
            std::cerr << "ERROR: invalid edge " << path[i - 1] << " " << path[i] << std::endl;
            return false;
        }

        seq += next;
    }

    assert(seq.size() == path.size() + graph.get_k() - 1 - offset);

    return true;
}

template <typename NodeType>
bool Alignment<NodeType>::is_valid(const DeBruijnGraph &graph,
                                   const DBGAlignerConfig *config) const {
    if (query_begin_ > query_end_) {
        std::cerr << "ERROR: query begin after query end" << std::endl
                  << *this << std::endl;
        return false;
    }

    std::string path;
    if (!spell_path(graph, nodes_, path, offset_)) {
        std::cerr << *this << std::endl;
        return false;
    }

    if (path != sequence_) {
        std::cerr << "ERROR: stored sequence is incorrect" << std::endl
                  << path << "\t"
                  << *this << std::endl;
        return false;
    }

    if (!cigar_.is_valid(sequence_, get_query())) {
        std::cerr << *this << std::endl;
        return false;
    }

    score_t cigar_score = config ? config->score_cigar(sequence_, get_query(), cigar_) : 0;
    if (config && score_ != cigar_score) {
        std::cerr << "ERROR: mismatch between CIGAR and score" << std::endl
                  << "CIGAR score: " << cigar_score << std::endl
                  << get_query() << "\t"
                  << *this << std::endl;
        return false;
    }

    return true;
}


template <typename NodeType>
QueryAlignment<NodeType>::QueryAlignment(const QueryAlignment &other)
      : query_(other.query_),
        query_rc_(other.query_rc_),
        alignments_(other.alignments_) {
    fix_pointers(other.get_query(),
                 other.get_query_reverse_complement());
}

template <typename NodeType>
QueryAlignment<NodeType>::QueryAlignment(QueryAlignment&& other) noexcept
      : query_(other.query_),
        query_rc_(other.query_rc_),
        alignments_(std::move(other.alignments_)) {
    fix_pointers(other.get_query(),
                 other.get_query_reverse_complement());
}

template <typename NodeType>
QueryAlignment<NodeType>::QueryAlignment(const std::string_view query) {
    // TODO: remove const_cast
    auto &qu = const_cast<std::string&>(query_);
    auto &qu_rc = const_cast<std::string&>(query_rc_);

    // pad sequences for easier access in 64-bit blocks
    qu.reserve(query.size() + 8);
    qu.resize(query.size());

    // TODO: use alphabet encoder
    // transform to upper and fix weird characters
    std::transform(query.begin(), query.end(), qu.begin(), [](char c) {
        return c >= 0 ? toupper(c) : 127;
    });
    memset(qu.data() + qu.size(), '\0', qu.capacity() - qu.size());

    qu_rc.reserve(query.size() + 8);
    qu_rc.resize(query.size());
    memcpy(qu_rc.data(), qu.data(), qu.capacity());
    reverse_complement(qu_rc.begin(), qu_rc.end());
}


template <typename NodeType>
void QueryAlignment<NodeType>
::fix_pointers(const std::string &query, const std::string &query_rc) {
    for (auto &alignment : alignments_) {
        const auto &other_query = alignment.get_orientation() ? query_rc : query;
        const auto &this_query = alignment.get_orientation() ? query_rc_ : query_;

        std::string_view query = alignment.get_query();
        assert(query.data() >= other_query.c_str());
        assert(query.data() + query.size() <= other_query.c_str() + other_query.size());

        alignment.set_query_begin(this_query.c_str() + (query.data() - other_query.c_str()));

        assert(alignment.get_query().begin() >= this_query.c_str());
        assert(alignment.get_query().end() <= this_query.c_str() + this_query.size());
    }
}


template class Alignment<>;
template class AlignmentPrefix<>;
template class AlignmentSuffix<>;
template class QueryAlignment<>;
template class DPTable<>;

} // namespace align
} // namespace graph
} // namespace mtg
