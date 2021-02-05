#include "aligner_alignment.hpp"

#include "aligner_dp_table.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;


template <typename NodeType>
Alignment<NodeType>::Alignment(std::string_view query,
                               std::vector<NodeType>&& nodes,
                               std::string&& sequence,
                               score_t score,
                               size_t clipping,
                               bool orientation,
                               size_t offset)
      : query_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
        score_(score), orientation_(orientation), offset_(offset) {
    size_t min_length = std::min(query_.size(), sequence_.size());

    cigar_ = std::inner_product(
        query_.data(),
        query_.data() + min_length,
        sequence_.c_str(),
        Cigar(Cigar::CLIPPED, clipping),
        [&](Cigar &cigar, bool equal) -> Cigar& {
            cigar.append(equal ? Cigar::MATCH : Cigar::MISMATCH);
            return cigar;
        },
        std::equal_to<char>()
    );

    assert(!(query_.size() - min_length) || (sequence_.size() - min_length));
    cigar_.append(Cigar::INSERTION, query_.size() - min_length);
    cigar_.append(Cigar::DELETION, sequence_.size() - min_length);
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const DPTable<NodeType> &dp_table,
                               const DBGAlignerConfig &config,
                               std::string_view query_view,
                               typename DPTable<NodeType>::const_iterator column,
                               size_t start_pos,
                               size_t offset,
                               NodeType *start_node,
                               const Alignment &seed)
      : orientation_(seed.get_orientation()), offset_(offset) {
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

    if (op == Cigar::DELETION)
        prev_node = prev_gap_node;

    uint32_t gap_count = op == Cigar::DELETION ? column->second.gap_count.at(i - shift) - 1 : 0;

    if (!i && prev_node == SequenceGraph::npos)
        return;

    // use config to recompute CIGAR score in DEBUG mode
    score_t score_track = score_;
    Cigar::Operator last_op = Cigar::CLIPPED;

    score_t gap_diff = config.gap_opening_penalty - config.gap_extension_penalty;

    // TODO: If there is a cyclic part of the graph in which the optimal
    //       alignment involves an insertion, then a score which was previously
    //       from a insertion may be replaced with a match. This will cause
    //       subsequent insertion scores to be wrong since they're no longer
    //       extensions. The only way to fix this is to store a separate vector
    //       to keep partial alignments ending in insertions.
    //       Until this is fixed, the score checking asserts have been commented out.

    std::vector<typename DPTable<NodeType>::const_iterator> out_columns;
    while (prev_node != SequenceGraph::npos) {
        auto prev_column = dp_table.find(prev_node);
        assert(prev_column != dp_table.end());
        assert(i || op == Cigar::DELETION);

        switch (op) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                --i;
                out_columns.emplace_back(column);

                if (last_op == Cigar::DELETION)
                    score_track -= gap_diff;

                // assert(column->second.scores.at(i + 1) >= score_track);
                score_track -= config.get_row(column->second.last_char)[query_view[i]];
                // assert(prev_column->second.scores.at(i) >= score_track);

            } break;
            case Cigar::DELETION: {
                out_columns.emplace_back(column);

                score_track -= config.gap_extension_penalty;

            } break;
            case Cigar::INSERTION: {
                assert(column == prev_column);
                --i;

                assert(i >= shift);
                assert(i - shift < column->second.scores.size());

                // assert(column->second.prev_nodes.at(i + 1) == 0xFF);
                // assert(column->second.scores.at(i + 1) >= score_track);

                if (column->second.ops.at(i - shift) == Cigar::DELETION) {
                    logger->error("INSERTION after DELETION: {}", query_view);
                    exit(1);
                }

                score_track -= column->second.ops.at(i - shift) == Cigar::INSERTION
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

            if (op == Cigar::DELETION)
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

        if (op == Cigar::DELETION)
            prev_node = prev_gap_node;
    }

    const auto &score_col = column->second.scores;

    if (last_op == Cigar::DELETION)
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

    query_ = { query_view.data() + i, start_pos - i };

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
        logger->warn(
            "Correcting score: {} -> {}\nQuery: {}\nTarget: {}\nSeed: {}\nExtension: {}",
            score_,
            score_ + correction,
            seed.get_sequence() + std::string(query_view),
            seed.get_sequence() + sequence_,
            seed.get_cigar().to_string(),
            cigar_.to_string()
        );
    }
}

template <typename NodeType>
void Alignment<NodeType>::append(Alignment&& other) {
    assert(query_.data() + query_.size() == other.query_.data());
    assert(orientation_ == other.orientation_);
    assert(cigar_.empty() || cigar_.back().first != Cigar::CLIPPED);

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;

    cigar_.append(std::move(other.cigar_));
    query_ = { query_.data(), query_.size() + other.query_.size() };
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
                || it->first == Cigar::INSERTION) {
            ++it;
            counter = 0;
            continue;
        }

        size_t jump = std::min({ offset_, static_cast<size_t>(it->second),
                                          static_cast<size_t>(nodes_.end() - jt) });
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
                                             std::string_view query_rev_comp) {
    assert(graph.is_canonical_mode());

    if (empty())
        return;

    assert(query_.size() + get_end_clipping() == query_rev_comp.size() - get_clipping());
    assert(is_valid(graph));

    trim_offset();
    assert(is_valid(graph));

    if (!offset_) {
        reverse_complement_seq_path(graph, sequence_, nodes_);
    } else {
        assert(nodes_.size() == 1);
        // extract target sequence prefix
        sequence_ = graph.get_node_sequence(nodes_[0]).substr(0, offset_) + sequence_;

        if (sequence_[0] == boss::BOSS::kSentinel) {
            // If the alignment starts from a source k-mer, then traverse forwards
            // until a non-dummy k-mer is hit and check if its reverse complement exists.

            // An appropriate node may not exist if offset_ is greater than the
            // number of sentinel characters (i.e., if some non-sentinel characters
            // from the node prefix are not included).

            const auto &dbg_succ = dynamic_cast<const DBGSuccinct&>(graph);
            const auto &boss = dbg_succ.get_boss();

#ifndef NDEBUG
            size_t num_sentinels = sequence_.find_last_of(boss::BOSS::kSentinel) + 1;
            assert(offset_ >= num_sentinels);
#endif

            // the node is present in the underlying graph, so use
            // lower-level methods
            boss::BOSS::edge_index edge = dbg_succ.kmer_to_boss_index(nodes_[0]);
            boss::BOSS::TAlphabet edge_label = boss.get_W(edge) % boss.alph_size;

            // TODO: This picks the node which is found by always traversing
            // the last outgoing edge. Is there a better way to pick a node?
            for (size_t i = 0; i < offset_; ++i) {
                edge = boss.fwd(edge, edge_label);
                edge_label = boss.get_W(edge) % boss.alph_size;
                if (edge_label == boss::BOSS::kSentinelCode) {
                    // reverse complement not found
                    assert(offset_ > num_sentinels);
                    *this = Alignment();
                    return;
                }

                nodes_[0] = dbg_succ.boss_to_kmer_index(edge);
                assert(nodes_[0]);
                sequence_.push_back(boss.decode(edge_label));
                assert(graph.get_node_sequence(nodes_[0])
                    == sequence_.substr(sequence_.size() - graph.get_k()));
            }

            assert(sequence_.size() == dbg_succ.get_k() + offset_);
            sequence_ = sequence_.substr(offset_);

            assert(nodes_ == map_sequence_to_nodes(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            sequence_.assign(sequence_.data() + offset_, graph.get_k() - offset_);

        } else {
            assert(nodes_.size() == 1);
            assert(nodes_ == map_sequence_to_nodes(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            // trim off ending from reverse complement (corresponding to the added prefix)
            for (size_t i = 0; i < offset_; ++i) {
                size_t indegree = 0;
                graph.adjacent_incoming_nodes(nodes_[0], [&](NodeType prev) {
                    ++indegree;

                    // TODO: there are multiple possible reverse complements, which
                    // do we pick? Currently we pick the first one
                    if (indegree == 1)
                        nodes_[0] = prev;
                });

                if (!indegree) {
                    *this = Alignment();
                    return;
                }

                sequence_.pop_back();
                assert(graph.get_node_sequence(nodes_[0]).substr(i + 1)
                    == sequence_.substr(sequence_.size() - graph.get_k() + i + 1));
            }

            assert(sequence_.size() == graph.get_k() - offset_);
        }

        assert(graph.get_node_sequence(nodes_[0]).substr(offset_) == sequence_);
    }

    std::reverse(cigar_.begin(), cigar_.end());
    assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

    orientation_ = !orientation_;
    query_ = { query_rev_comp.data() + get_clipping(),
               query_rev_comp.size() - get_clipping() - get_end_clipping() };
    assert(is_valid(graph));
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
template <typename NodeType>
Json::Value Alignment<NodeType>::path_json(size_t node_size,
                                           std::string_view label) const {
    assert(nodes_.size());

    Json::Value path;

    auto cigar_it = cigar_.begin();
    if (cigar_.size() && cigar_it->first == Cigar::CLIPPED) {
        cigar_it++;
    }

    size_t cigar_offset = 0;
    assert(cigar_it != cigar_.end());

    int64_t rank = 1;
    const char *query_start = query_.data();

#ifndef NDEBUG
    const char *query_end = query_.data() + query_.size();
#endif

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
                assert(query_start + next_size <= query_end);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::INSERTION: {
                assert(query_start + next_size <= query_end);
                // this assumes that INSERTIONs can't happen right after DELETIONs
                //edit["from_length"] = 0;
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;

                // the target is not consumed, so reset the position
                next_pos = cur_pos;
            } break;
            case Cigar::DELETION: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                assert(query_start + next_size <= query_end);
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

        if (cigar_it->first == Cigar::INSERTION) {
            Json::Value edit;
            size_t length = cigar_it->second - cigar_offset;
            assert(query_start + length < query_end);
            // TODO: this assumes that INSERTIONs can't happen right after DELETIONs
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
                assert(query_start < query_end);
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                edit["sequence"] = std::string(query_start, 1);
                query_start++;
            } break;
            case Cigar::DELETION: {
                edit["from_length"] = 1;
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                query_start++;
            } break;
            case Cigar::INSERTION:
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

    assert(query_start == query_end);
    assert(cigar_it == cigar_.end()
            || (cigar_it + 1 == cigar_.end() && cigar_it->first == Cigar::CLIPPED));

    path["length"] = Json::Value::UInt64(nodes_.size());

    if (nodes_.front() == nodes_.back())
        path["is_circular"] = true;

    if (label.data())
        path["name"] = std::string(label);

    return path;
}

template <typename NodeType>
Json::Value Alignment<NodeType>::to_json(std::string_view full_query,
                                         const DeBruijnGraph &graph,
                                         bool is_secondary,
                                         std::string_view read_name,
                                         std::string_view label) const {
    assert(is_valid(graph));

    // encode alignment
    Json::Value alignment;

    alignment["name"] = read_name.data() ? std::string(read_name) : "";
    alignment["sequence"] = std::string(full_query);

    if (sequence_.size())
        alignment["annotation"]["ref_sequence"] = sequence_;

    if (query_.empty())
        return alignment;

    assert(query_.data() - cigar_.get_clipping() == full_query.data());
    assert(query_.size() + cigar_.get_clipping() + cigar_.get_end_clipping()
            == full_query.size());

    alignment["annotation"]["cigar"] = cigar_.to_string();

    // encode path
    if (nodes_.size())
        alignment["path"] = path_json(graph.get_k(), label);

    alignment["score"] = static_cast<int32_t>(score_);

    if (cigar_.get_clipping()) {
        alignment["query_position"] = static_cast<int32_t>(cigar_.get_clipping());
        alignment["soft_clipped"] = static_cast<bool>(cigar_.get_clipping());
    }

    if (is_secondary)
        alignment["is_secondary"] = is_secondary;

    alignment["identity"] = query_.size()
        ? static_cast<double>(get_num_matches()) / query_.size()
        : 0;

    alignment["read_mapped"] = static_cast<bool>(query_.size());

    if (orientation_)
        alignment["read_on_reverse_strand"] = orientation_;

    if (label.data())
        alignment["sample_name"] = std::string(label);


    // Unused flags (for now)
    //alignment["quality"]; // bytes
    //alignment["mapping_quality"]; // int32
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
::load_from_json(const Json::Value &alignment, const DeBruijnGraph &graph) {
    cigar_.clear();
    nodes_.clear();
    sequence_.clear();

    auto query_sequence = std::make_shared<const std::string>(
        alignment["sequence"].asString()
    );
    orientation_ = alignment["read_on_reverse_strand"].asBool();
    score_ = alignment["score"].asInt();
    const Json::Value &mapping = alignment["path"]["mapping"];
    assert(mapping.size() == alignment["path"]["length"].asUInt64());

    cigar_.append(Cigar::CLIPPED, alignment["query_position"].asInt());

    const char *this_query_begin = query_sequence->c_str() + get_clipping();
    size_t this_query_size = 0;

    size_t path_steps = 0;
    Json::ArrayIndex i = 0;
    offset_ = mapping[i]["position"]["offset"].asUInt64();

    for ( ; i < mapping.size(); ++i) {
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
        const Json::Value &edits = mapping[i]["edit"];

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
                this_query_size += edits[j]["to_length"].asUInt64();

            } else if (edits[j]["from_length"].asUInt64()) {
                cigar_.append(Cigar::DELETION, edits[j]["from_length"].asUInt64());
                path_steps += edits[j]["from_length"].asUInt64();

            } else {
                assert(edits[j]["to_length"].asUInt64());
                cigar_.append(Cigar::INSERTION, edits[j]["to_length"].asUInt64());
                this_query_size += edits[j]["to_length"].asUInt64();
            }
        }
    }

    query_ = { this_query_begin, this_query_size };

    if (query_.size() + get_clipping() < query_sequence->size()) {
        cigar_.append(Cigar::CLIPPED,
                      query_sequence->size() - query_.size() - get_clipping());
    }

    if (alignment["annotation"]["cigar"]
            && cigar_ != Cigar(alignment["annotation"]["cigar"].asString())) {
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

    if (!cigar_.is_valid(sequence_, query_)) {
        std::cerr << *this << std::endl;
        return false;
    }

    score_t cigar_score = config ? config->score_cigar(sequence_, query_, cigar_) : 0;
    if (config && score_ != cigar_score) {
        std::cerr << "ERROR: mismatch between CIGAR and score" << std::endl
                  << "CIGAR score: " << cigar_score << std::endl
                  << query_ << "\t"
                  << *this << std::endl;
        return false;
    }

    return true;
}


template <typename NodeType>
QueryAlignment<NodeType>::QueryAlignment(std::string_view query,
                                         bool is_reverse_complement)
          : query_(new std::string()), query_rc_(new std::string()) {
    // pad sequences for easier access in 64-bit blocks
    query_->reserve(query.size() + 8);
    query_->resize(query.size());
    query_rc_->reserve(query.size() + 8);
    query_rc_->resize(query.size());

    // TODO: use alphabet encoder
    // transform to upper and fix non-standard characters
    std::transform(query.begin(), query.end(), query_->begin(), [](char c) {
        return c >= 0 ? toupper(c) : 127;
    });

    // fill padding with '\0'
    memset(query_->data() + query_->size(), '\0', query_->capacity() - query_->size());

    // set the reverse complement
    memcpy(query_rc_->data(), query_->data(), query_->capacity());
    reverse_complement(query_rc_->begin(), query_rc_->end());

    if (is_reverse_complement)
        std::swap(query_, query_rc_);
}


template class Alignment<>;
template struct LocalAlignmentLess<>;
template class QueryAlignment<>;

} // namespace align
} // namespace graph
} // namespace mtg
