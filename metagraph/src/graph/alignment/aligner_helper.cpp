#include "aligner_helper.hpp"

#include <tsl/hopscotch_set.h>

#include "common/algorithms.hpp"
#include "common/logger.hpp"

using mg::common::logger;


template <typename NodeType>
bool DPTable<NodeType>::add_seed(NodeType start_node,
                                 char start_char,
                                 score_t initial_score,
                                 score_t min_score,
                                 size_t size,
                                 size_t start_pos,
                                 int8_t gap_opening_penalty,
                                 int8_t gap_extension_penalty,
                                 size_t query_offset) {
    assert(start_pos < size);

    query_offset_ = query_offset;
    start_node_ = start_node;

    // Initialize first column
    auto &table_init = dp_table_[start_node];
    if (!table_init.size())
        table_init = Column(size, min_score, start_char, start_pos);

    bool update = false;

    if (table_init.best_score() < initial_score) {
        table_init.scores[start_pos] = initial_score;
        table_init.ops[start_pos] = Cigar::Operator::MATCH;
        table_init.prev_nodes[start_pos] = SequenceGraph::npos;
        update = true;
    }

    if (size - start_pos <= 1)
        return update;

    if (table_init.scores[start_pos] + gap_opening_penalty > table_init.scores[start_pos + 1]) {
        table_init.scores[start_pos + 1] = table_init.scores[start_pos] + gap_opening_penalty;
        table_init.ops[start_pos + 1] = Cigar::Operator::DELETION;
        table_init.prev_nodes[start_pos + 1] = start_node;
        update = true;
    }

    for (size_t i = start_pos + 2; i < size; ++i) {
        if (table_init.scores[i - 1] + gap_extension_penalty <= table_init.scores[i])
            break;

        table_init.scores[i] = table_init.scores[i - 1] + gap_extension_penalty;
        table_init.ops[i] = Cigar::Operator::DELETION;
        table_init.prev_nodes[i] = start_node;
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
                     bool orientation,
                     score_t min_path_score,
                     NodeType *node) {
    NodeType start_node;
    if (config.num_alternative_paths == 1 && node) {
        // avoid sorting column iterators if we're only interested in the top path_
        start_node = DeBruijnGraph::npos;
        auto column_it = dp_table_.find(*node);
        assert(column_it != dp_table_.end());
        Alignment<NodeType> alignment(*this,
                                      query_view,
                                      column_it,
                                      column_it->second.best_pos,
                                      orientation,
                                      graph.get_k() - 1,
                                      &start_node);

        if (UNLIKELY(alignment.empty() && !alignment.get_query().data())) {
            return;
        }

        // TODO: remove this when the branch and bound is set to only consider
        //       converged scores
        alignment.recompute_score(config);

        assert(alignment.is_valid(graph, &config));

        callback(std::move(alignment), start_node);
    }

    // store visited nodes in paths to avoid returning subalignments
    tsl::hopscotch_set<key_type> visited_nodes;

    std::vector<const_iterator> starts;
    starts.reserve(size());
    for (const_iterator it = dp_table_.begin(); it != dp_table_.end(); ++it) {
        if (it->second.best_score() > min_path_score
                && it->second.best_op() == Cigar::Operator::MATCH)
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
                                 query_view,
                                 column_it,
                                 column_it->second.best_pos,
                                 orientation,
                                 graph.get_k() - 1,
                                 &start_node);

        if (UNLIKELY(next.empty() && !next.get_query().data())) {
            continue;
        }

        // TODO: remove this when the branch and bound is set to only consider
        //       converged scores
        next.recompute_score(config);
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
        Cigar(Cigar::Operator::CLIPPED, clipping),
        [&](Cigar &cigar, bool equal) -> Cigar& {
            cigar.append(equal
                  ? Cigar::Operator::MATCH
                  : Cigar::Operator::MISMATCH
              );
            return cigar;
        },
        std::equal_to<char>()
    );

    assert(!(query_size - min_length) || (sequence_.size() - min_length));
    cigar_.append(Cigar::Operator::DELETION, query_size - min_length);
    cigar_.append(Cigar::Operator::INSERTION, sequence_.size() - min_length);
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const DPTable &dp_table,
                               const std::string_view query_view,
                               typename DPTable::const_iterator column,
                               size_t start_pos,
                               bool orientation,
                               size_t offset,
                               NodeType *start_node)
      : query_begin_(NULL),
        query_end_(NULL),
        score_(column->second.scores.at(start_pos)),
        orientation_(orientation),
        offset_(offset) {
    assert(start_node);

    auto i = start_pos;
    const auto* op = &column->second.ops.at(i);
    const auto* prev_node = &column->second.prev_nodes.at(i);

    if (!i && *prev_node == SequenceGraph::npos)
        return;

    std::vector<typename DPTable::const_iterator> out_columns;
    while (*prev_node != SequenceGraph::npos) {
        cigar_.append(*op);

        if (*op != Cigar::Operator::DELETION)
            out_columns.emplace_back(column);

        if (*op != Cigar::Operator::INSERTION)
            --i;

        column = dp_table.find(*prev_node);
        op = &column->second.ops.at(i);
        prev_node = &column->second.prev_nodes.at(i);
    }

    if (column->first == dp_table.get_start_node())
        score_ -= column->second.scores.front();

    *start_node = column->first;

    if (UNLIKELY(i > std::numeric_limits<Cigar::LengthType>::max())) {
        throw std::runtime_error("Error: clipping length can't be stored in CIGAR");
    }

    cigar_.append(Cigar::Operator::CLIPPED, i);
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
}

template <typename NodeType>
void Alignment<NodeType>::append(Alignment&& other) {
    assert(query_end_ == other.query_begin_);
    assert(orientation_ == other.orientation_);
    assert(cigar_.empty() || cigar_.back().first != Cigar::Operator::CLIPPED);

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;

    cigar_.append(std::move(other.cigar_));
    query_end_ = other.query_end_;
}

template <typename NodeType>
void Alignment<NodeType>::recompute_score(const DBGAlignerConfig &config) {
    auto new_score = config.score_cigar(sequence_, get_query(), cigar_);

    if (utils::get_verbose() && score_ != new_score)
        logger->trace("changing score from {}, to {}", score_, new_score);

    score_ = new_score;
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
template <typename NodeType>
Json::Value Alignment<NodeType>::path_json(size_t node_size,
                                           const std::string &label) const {
    assert(nodes_.size());

    Json::Value path;

    auto cigar_it = cigar_.begin();
    if (cigar_.size() && cigar_it->first == Cigar::Operator::CLIPPED) {
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

    while (cur_pos < node_size && cigar_it != cigar_.end()) {
        assert(cigar_it->second > cigar_offset);
        size_t next_pos = std::min(node_size,
                                   cur_pos + (cigar_it->second - cigar_offset));
        size_t next_size = next_pos - cur_pos;
        assert(cigar_offset + next_size <= cigar_it->second);

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::Operator::MISMATCH: {
                assert(query_start + next_size <= query_end_);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::DELETION: {
                assert(query_start + next_size <= query_end_);
                // this assumes that DELETIONS can't happen right after insertions
                //edit["from_length"] = 0;
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::INSERTION: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                //edit["to_length"] = 0;
            } break;
            case Cigar::Operator::MATCH: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::CLIPPED: {
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

        if (cigar_it->first == Cigar::Operator::DELETION) {
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
            case Cigar::Operator::MISMATCH: {
                assert(query_start < query_end_);
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                edit["sequence"] = std::string(query_start, 1);
                query_start++;
            } break;
            case Cigar::Operator::INSERTION: {
                edit["from_length"] = 1;
                //edit["to_length"] = 0;
            } break;
            case Cigar::Operator::MATCH: {
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                query_start++;
            } break;
            case Cigar::Operator::DELETION:
            case Cigar::Operator::CLIPPED: assert(false); break;
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
            || (cigar_it + 1 == cigar_.end()
                    && cigar_it->first == Cigar::Operator::CLIPPED));

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

    cigar_.append(Cigar::Operator::CLIPPED, query_begin_ - query_start);
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
                    cigar_.append(Cigar::Operator::MISMATCH,
                                  edits[j]["from_length"].asUInt64());
                } else {
                    cigar_.append(Cigar::Operator::MATCH,
                                  edits[j]["from_length"].asUInt64());
                }

                path_steps += edits[j]["from_length"].asUInt64();
                query_end_ += edits[j]["to_length"].asUInt64();
            } else if (edits[j]["from_length"].asUInt64()) {
                cigar_.append(Cigar::Operator::INSERTION,
                              edits[j]["from_length"].asUInt64());

                path_steps += edits[j]["from_length"].asUInt64();
            } else {
                assert(edits[j]["to_length"].asUInt64());
                cigar_.append(Cigar::Operator::DELETION,
                              edits[j]["to_length"].asUInt64());

                query_end_ += edits[j]["to_length"].asUInt64();
            }
        }
    }

    if (query_end_ < query_sequence->c_str() + query_sequence->length()) {
        cigar_.append(Cigar::Operator::CLIPPED,
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
                  << *this << std::endl;
        return false;
    }

    if (config && score_ != config->score_cigar(sequence_, get_query(), cigar_)) {
        std::cerr << "ERROR: mismatch between CIGAR and score" << std::endl
                  << *this << std::endl;
        return false;
    }

    if (!cigar_.is_valid(sequence_, get_query())) {
        std::cerr << *this << std::endl;
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
QueryAlignment<NodeType>::QueryAlignment(const std::string &query)
      : query_(query),
        query_rc_(query) {
    // TODO: remove const_cast
    reverse_complement(const_cast<char*>(query_rc_.data()),
                       const_cast<char*>(query_rc_.data() + query_rc_.size()));
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
template class QueryAlignment<>;
template class DPTable<>;
