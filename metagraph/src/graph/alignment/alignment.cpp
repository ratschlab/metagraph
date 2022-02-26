#include "alignment.hpp"

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/rc_dbg.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

Alignment::Alignment(std::string_view query,
                     std::vector<node_index>&& nodes,
                     std::string&& sequence,
                     score_t score,
                     size_t clipping,
                     bool orientation,
                     size_t offset)
      : query_view_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
        score_(score),
        orientation_(orientation),
        offset_(offset) {
    assert(query_view_.size() == sequence_.size());
    cigar_ = std::inner_product(query_view_.begin(), query_view_.end(), sequence_.begin(),
                                Cigar(Cigar::CLIPPED, clipping),
                                [&](Cigar &cigar, bool equal) -> Cigar& {
                                    cigar.append(equal ? Cigar::MATCH : Cigar::MISMATCH);
                                    return cigar;
                                },
                                std::equal_to<char>());
}

std::string Alignment::format_coords() const {
    if (!label_coordinates.size())
        return "";

    assert(label_columns.size());
    assert(label_encoder);
    assert(label_coordinates.size() == label_columns.size());

    std::vector<std::string> decoded_labels;
    decoded_labels.reserve(label_columns.size());

    for (size_t i = 0; i < label_columns.size(); ++i) {
        decoded_labels.emplace_back(label_encoder->decode(label_columns[i]));
        for (uint64_t coord : label_coordinates[i]) {
            // alignment coordinates are 1-based inclusive ranges
            decoded_labels.back()
                += fmt::format(":{}-{}", coord + 1, coord + sequence_.size());
        }
    }

    return fmt::format("{}", fmt::join(decoded_labels, ";"));
}

bool Alignment::append(Alignment&& other) {
    assert(query_view_.data() + query_view_.size() + other.get_clipping()
            == other.query_view_.data());
    assert(orientation_ == other.orientation_);

    bool ret_val = false;

    if (label_coordinates.size() && other.label_coordinates.empty())
        label_coordinates.clear();

    if (label_columns.size() && other.label_columns.empty())
        label_columns.clear();

    if (label_coordinates.size()) {
        assert(label_columns.size() == label_coordinates.size());
        Columns merged_label_columns;
        CoordinateSet merged_label_coordinates;

        // if the alignments fit together without gaps, make sure that the
        // coordinates form a contiguous range
        CoordIntersection coord_inter(sequence_.size());
        utils::match_indexed_values(
            label_columns.begin(), label_columns.end(),
            label_coordinates.begin(),
            other.label_columns.begin(), other.label_columns.end(),
            other.label_coordinates.begin(),
            [&](auto col, const auto &coords, const auto &other_coords) {
                Tuple merged;
                coord_inter(coords.begin(), coords.end(),
                            other_coords.begin(), other_coords.end(),
                            std::back_inserter(merged));
                if (merged.size()) {
                    merged_label_columns.push_back(col);
                    merged_label_coordinates.push_back(std::move(merged));
                }
            }
        );

        if (merged_label_columns.empty()) {
            *this = Alignment();
            return true;
        }

        ret_val = merged_label_columns.size() < label_columns.size();

        if (!ret_val) {
            for (size_t i = 0; i < label_columns.size(); ++i) {
                if (merged_label_coordinates[i].size() < label_coordinates[i].size()) {
                    ret_val = true;
                    break;
                }
            }
        }

        std::swap(label_columns, merged_label_columns);
        std::swap(label_coordinates, merged_label_coordinates);

    } else if (label_columns.size()) {
        Columns merged_label_columns;
        std::set_intersection(label_columns.begin(), label_columns.end(),
                              other.label_columns.begin(), other.label_columns.end(),
                              std::back_inserter(merged_label_columns));

        if (merged_label_columns.empty()) {
            *this = Alignment();
            return true;
        }

        ret_val = merged_label_columns.size() < label_columns.size();

        std::swap(label_columns, merged_label_columns);
    }

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;
    cigar_.append(std::move(other.cigar_));
    // expand the query window to cover both alignments
    query_view_ = std::string_view(query_view_.data(),
                                   other.query_view_.end() - query_view_.begin());
    return ret_val;
}

size_t Alignment::trim_offset() {
    if (!offset_ || nodes_.size() <= 1)
        return 0;

    assert(nodes_.front());

    size_t first_dummy = (std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
        - nodes_.begin()) - 1;
    size_t trim = std::min(std::min(offset_, nodes_.size() - 1), first_dummy);
    offset_ -= trim;
    nodes_.erase(nodes_.begin(), nodes_.begin() + trim);
    assert(nodes_.front());
    return trim;
}

size_t Alignment::trim_query_prefix(size_t n,
                                    size_t node_overlap,
                                    const DBGAlignerConfig &config,
                                    bool trim_excess_deletions) {
    size_t clipping = get_clipping();
    const char *query_begin = query_view_.data() - clipping;
    auto it = cigar_.data().begin() + static_cast<bool>(clipping);
    size_t cigar_offset = 0;

    auto s_it = sequence_.begin();
    auto node_it = nodes_.begin();

    auto consume_ref = [&]() {
        assert(s_it != sequence_.end());
        ++s_it;
        if (offset_ < node_overlap) {
            ++offset_;
        } else if (node_it + 1 < nodes_.end()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_deletions && it->first == Cigar::DELETION)) {
        if (it == cigar_.data().end()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.end());
                score_ -= config.score_matrix[query_view_[0]][*s_it];
                query_view_.remove_prefix(1);
                --n;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_prefix(1);
                --n;
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    size_t seq_trim = s_it - sequence_.begin();
    for (auto &coords : label_coordinates) {
        for (auto &c : coords) {
            c += seq_trim;
        }
    }

    if (!clipping && it != cigar_.data().begin())
        score_ -= config.left_end_bonus;

    nodes_.erase(nodes_.begin(), node_it);
    sequence_.erase(sequence_.begin(), s_it);
    it->second -= cigar_offset;
    cigar_.data().erase(cigar_.data().begin(), it);
    extend_query_begin(query_begin);

    return cigar_offset;
}

size_t Alignment::trim_query_suffix(size_t n,
                                    const DBGAlignerConfig &config,
                                    bool trim_excess_deletions) {
    size_t end_clipping = get_end_clipping();
    const char *query_end = query_view_.data() + query_view_.size() + end_clipping;

    trim_end_clipping();
    auto it = cigar_.data().rbegin();
    size_t cigar_offset = 0;

    auto s_it = sequence_.rbegin();
    auto node_it = nodes_.rbegin();

    auto consume_ref = [&]() {
        assert(s_it != sequence_.rend());
        ++s_it;
        if (node_it + 1 < nodes_.rend()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_deletions && it->first == Cigar::DELETION)) {
        if (it == cigar_.data().rend()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.rend());
                score_ -= config.score_matrix[query_view_.back()][*s_it];
                query_view_.remove_suffix(1);
                --n;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_suffix(1);
                --n;
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    if (!end_clipping && (cigar_offset || it.base() != cigar_.data().end()))
        score_ -= config.right_end_bonus;

    nodes_.erase(node_it.base(), nodes_.end());
    sequence_.erase(s_it.base(), sequence_.end());
    it->second -= cigar_offset;
    cigar_.data().erase(it.base(), cigar_.data().end());

    extend_query_end(query_end);

    return cigar_offset;
}

size_t Alignment::trim_reference_prefix(size_t n,
                                        size_t node_overlap,
                                        const DBGAlignerConfig &config,
                                        bool trim_excess_insertions) {
    size_t clipping = get_clipping();
    const char *query_begin = query_view_.data() - clipping;

    auto it = cigar_.data().begin() + static_cast<bool>(clipping);
    size_t cigar_offset = 0;

    auto s_it = sequence_.begin();
    auto node_it = nodes_.begin();

    size_t seq_trim = 0;

    auto consume_ref = [&]() {
        assert(s_it != sequence_.end());
        assert(n);
        if (*s_it != '$') {
            ++seq_trim;
            --n;
        }

        ++s_it;
        if (offset_ < node_overlap) {
            ++offset_;
        } else if (node_it + 1 < nodes_.end()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_insertions && it->first == Cigar::INSERTION)) {
        if (it == cigar_.data().end()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.end());
                score_ -= config.score_matrix[query_view_[0]][*s_it];
                query_view_.remove_prefix(1);
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_prefix(1);
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::NODE_INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
            } break;
            case Cigar::CLIPPED: {
                assert(false && "this should not happen");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    for (auto &coords : label_coordinates) {
        for (auto &c : coords) {
            c += seq_trim;
        }
    }

    if (!clipping && (cigar_offset || it != cigar_.data().begin()))
        score_ -= config.left_end_bonus;

    nodes_.erase(nodes_.begin(), node_it);
    sequence_.erase(sequence_.begin(), s_it);
    it->second -= cigar_offset;
    cigar_.data().erase(cigar_.data().begin(), it);

    extend_query_begin(query_begin);

    return cigar_offset;
}

size_t Alignment::trim_reference_suffix(size_t n,
                                        const DBGAlignerConfig &config,
                                        bool trim_excess_insertions) {
    size_t end_clipping = get_end_clipping();
    const char *query_end = query_view_.data() + query_view_.size() + end_clipping;

    trim_end_clipping();
    auto it = cigar_.data().rbegin();
    size_t cigar_offset = 0;

    auto s_it = sequence_.rbegin();
    auto node_it = nodes_.rbegin();

    auto consume_ref = [&]() {
        --n;
        assert(s_it != sequence_.rend());
        ++s_it;
        if (node_it + 1 < nodes_.rend()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_insertions && it->first == Cigar::INSERTION)) {
        if (it == cigar_.data().rend()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.rend());
                score_ -= config.score_matrix[query_view_.back()][*s_it];
                query_view_.remove_suffix(1);
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_suffix(1);
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    if (!end_clipping && (cigar_offset || it.base() != cigar_.data().end()))
        score_ -= config.right_end_bonus;

    nodes_.erase(node_it.base(), nodes_.end());
    sequence_.erase(s_it.base(), sequence_.end());
    it->second -= cigar_offset;
    cigar_.data().erase(it.base(), cigar_.data().end());

    extend_query_end(query_end);

    return cigar_offset;
}

void Alignment::reverse_complement(const DeBruijnGraph &graph,
                                   std::string_view query_rev_comp) {
    assert(query_view_.size() + get_end_clipping() == query_rev_comp.size() - get_clipping());

    trim_offset();
    assert(!offset_ || nodes_.size() == 1);

    if (dynamic_cast<const RCDBG*>(&graph)) {
        if (offset_) {
            *this = Alignment();
        } else {
            std::reverse(cigar_.data().begin(), cigar_.data().end());
            std::reverse(nodes_.begin(), nodes_.end());
            ::reverse_complement(sequence_.begin(), sequence_.end());
            assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

            orientation_ = !orientation_;
            query_view_ = { query_rev_comp.data() + get_clipping(),
                            query_rev_comp.size() - get_clipping() - get_end_clipping() };
        }
        return;
    }

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

            // TODO: this cascade of graph unwrapping is ugly, find a cleaner way to do it
            const DeBruijnGraph *base_graph = &graph;
            if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(base_graph))
                base_graph = &rc_dbg->get_graph();

            const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_graph);
            if (canonical)
                base_graph = &canonical->get_graph();

            const auto &dbg_succ = dynamic_cast<const DBGSuccinct&>(*base_graph);

            size_t num_sentinels = sequence_.find_last_of(boss::BOSS::kSentinel) + 1;
            assert(offset_ >= num_sentinels);

            if (canonical && nodes_[0] != canonical->get_base_node(nodes_[0])) {
                // reverse complement of a sink dummy k-mer, no point in traversing
                assert(num_sentinels == 1);
                *this = Alignment();
                return;
            }

            size_t num_first_steps = canonical ? std::min(offset_, num_sentinels) : offset_;

            // the node is present in the underlying graph, so use
            // lower-level methods
            const auto &boss = dbg_succ.get_boss();
            boss::BOSS::edge_index edge = dbg_succ.kmer_to_boss_index(nodes_[0]);
            boss::BOSS::TAlphabet edge_label = boss.get_W(edge) % boss.alph_size;

            // TODO: This picks the node which is found by always traversing
            // the last outgoing edge. Is there a better way to pick a node?
            for (size_t i = 0; i < num_first_steps; ++i) {
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

            for (size_t i = num_first_steps; i < offset_; ++i) {
                node_index next_node = 0;
                char last_char;
                canonical->call_outgoing_kmers(nodes_[0], [&](node_index next, char c) {
                    if (c == boss::BOSS::kSentinel)
                        return;

                    next_node = next;
                    last_char = c;
                });

                if (!next_node) {
                    *this = Alignment();
                    return;
                } else {
                    nodes_[0] = next_node;
                    sequence_.push_back(last_char);
                    assert(graph.get_node_sequence(nodes_[0])
                        == sequence_.substr(sequence_.size() - graph.get_k()));
                }
            }

            assert(sequence_.size() == dbg_succ.get_k() + offset_);
            sequence_ = sequence_.substr(offset_);

            assert(nodes_ == map_to_nodes_sequentially(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            sequence_.assign(sequence_.data() + offset_, graph.get_k() - offset_);

        } else {
            assert(nodes_.size() == 1);
            assert(nodes_ == map_to_nodes_sequentially(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            // trim off ending from reverse complement (corresponding to the added prefix)
            for (size_t i = 0; i < offset_; ++i) {
                size_t indegree = 0;
                graph.adjacent_incoming_nodes(nodes_[0], [&](node_index prev) {
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

    std::reverse(cigar_.data().begin(), cigar_.data().end());
    assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

    orientation_ = !orientation_;
    query_view_ = { query_rev_comp.data() + get_clipping(),
                    query_rev_comp.size() - get_clipping() - get_end_clipping() };
    assert(is_valid(graph));
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
Json::Value path_json(const std::vector<DeBruijnGraph::node_index> &nodes,
                      const Cigar &cigar,
                      size_t node_size,
                      std::string_view query_view,
                      size_t offset,
                      const std::string &label) {
    assert(nodes.size());

    Json::Value path;

    auto cigar_it = cigar.data().begin();
    if (cigar.size() && cigar_it->first == Cigar::CLIPPED) {
        cigar_it++;
    }

    size_t cigar_offset = 0;
    assert(cigar_it != cigar.data().end());

    int64_t rank = 1;
    const char *query_start = query_view.data();

#ifndef NDEBUG
    const char *query_end = query_view.data() + query_view.size();
#endif

    size_t cur_pos = rank == 1 ? offset : 0;

    Json::Value mapping;
    Json::Value position;
    position["node_id"] = nodes.front();

    if (cur_pos)
        position["offset"] = Json::Value::UInt64(cur_pos);

    // set to true if the node is the reverse complement of the query
    //position["is_reverse"] = false;

    mapping["position"] = position;

    // handle alignment to the first node
    while (cur_pos < node_size && cigar_it != cigar.data().end()) {
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
                assert(cigar_it == cigar.data().end());
                continue;
            } break;
            case Cigar::NODE_INSERTION: {
                assert(false && "this should not be reached");
            } break;
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
    for (auto node_it = nodes.begin() + 1; node_it != nodes.end(); ++node_it) {
        assert(cigar_it != cigar.data().end());
        assert(cigar_it->second > cigar_offset);

        Json::Value mapping;
        Json::Value position;
        position["node_id"] = *node_it;
        position["offset"] = Json::Value::UInt64(node_size - 1);
        // set to true if the node is the reverse complement of the query
        //position["is_reverse"] = false;
        mapping["position"] = position;

        if (cigar_it->first == Cigar::INSERTION || cigar_it->first == Cigar::CLIPPED) {
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
            assert(cigar_it != cigar.data().end());
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
            case Cigar::NODE_INSERTION:
            case Cigar::CLIPPED: assert(false && "this should not be reached"); break;
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
    assert(cigar_it == cigar.data().end()
            || (cigar_it + 1 == cigar.data().end() && cigar_it->first == Cigar::CLIPPED));

    path["length"] = Json::Value::UInt64(nodes.size());

    if (nodes.front() == nodes.back())
        path["is_circular"] = true;

    path["name"] = label;

    return path;
}

Json::Value Alignment::to_json(size_t node_size,
                               bool is_secondary,
                               const std::string &read_name,
                               const std::string &label) const {
    if (sequence_.find("$") != std::string::npos
            || std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos) != nodes_.end()) {
        throw std::runtime_error("JSON output for chains not supported");
    }

    std::string_view full_query = {
        query_view_.data() - get_clipping(),
        query_view_.size() + get_clipping() + get_end_clipping()
    };

    // encode alignment
    Json::Value alignment;

    alignment["name"] = read_name;
    alignment["sequence"] = std::string(full_query);

    if (sequence_.size())
        alignment["annotation"]["ref_sequence"] = sequence_;

    if (query_view_.empty())
        return alignment;

    assert(query_view_.data() - cigar_.get_clipping() == full_query.data());
    assert(query_view_.size() + cigar_.get_clipping() + cigar_.get_end_clipping()
            == full_query.size());

    alignment["annotation"]["cigar"] = cigar_.to_string();

    // encode path
    if (nodes_.size())
        alignment["path"] = path_json(nodes_, cigar_, node_size, query_view_, offset_, label);

    alignment["score"] = static_cast<int32_t>(score_);

    if (cigar_.get_clipping()) {
        alignment["query_position"] = static_cast<int32_t>(cigar_.get_clipping());
        alignment["soft_clipped"] = static_cast<bool>(cigar_.get_clipping());
    }

    if (is_secondary)
        alignment["is_secondary"] = is_secondary;

    alignment["identity"] = query_view_.size()
        ? static_cast<double>(cigar_.get_num_matches()) / query_view_.size()
        : 0;

    alignment["read_mapped"] = static_cast<bool>(query_view_.size());

    if (orientation_)
        alignment["read_on_reverse_strand"] = orientation_;

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

void Alignment::load_from_json(const Json::Value &alignment,
                               const DeBruijnGraph &graph,
                               std::string *query_sequence) {
    assert(query_sequence);

    cigar_ = Cigar();
    nodes_.clear();
    sequence_.clear();

    *query_sequence = alignment["sequence"].asString();
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
            graph.call_outgoing_kmers(*(nodes_.rbegin() + 1),
                                      [&](auto node, char c) {
                if (node == nodes_.back())
                    sequence_.push_back(c);
            });
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

    query_view_ = { this_query_begin, this_query_size };

    if (query_view_.size() + get_clipping() < query_sequence->size()) {
        cigar_.append(Cigar::CLIPPED,
                      query_sequence->size() - query_view_.size() - get_clipping());
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
}

void Alignment::insert_gap_prefix(ssize_t gap_length,
                                  size_t node_overlap,
                                  const DBGAlignerConfig &config) {
    size_t extra_nodes = node_overlap + 1;

    if (gap_length < 0) {
        // alignments overlap
        // extra_nodes = k - 1 - matching_overlap
        // e.g.,
        // k = 4
        // overlap = 4
        // matching overlap = 2
        // ATGCTATGCA
        //       ACCAACGACT

        trim_clipping();
        assert(extra_nodes + gap_length > 1);
        extra_nodes += gap_length - 1;

        if (offset_) {
            // if there are suffix-mapped nodes, only keep the ones that are
            // part of the overlap
            assert(static_cast<ssize_t>(offset_) >= -gap_length);
            nodes_.erase(nodes_.begin(), nodes_.begin() + offset_ + gap_length);
        }

        if (extra_nodes) {
            // they can't be joined since the overlap is too small
            // ATGCTATGCA
            //           ACGACT
            //       TGCA
            //        GCAA - added
            //         CAAC
            //          AACG
            //           ACGA
            score_ += config.gap_opening_penalty
                + (extra_nodes - 1) * config.gap_extension_penalty;
            cigar_.data().insert(cigar_.data().begin(),
                                 Cigar::value_type{ Cigar::NODE_INSERTION, extra_nodes });
        }
    } else {
        // no overlap
        // extra_nodes = k
        // e.g.,
        // k = 4
        // gap = 2
        // ATGCTATGCA
        //             ACGTACGACT
        //       TGCA
        //        GCA$ - added
        //         CA$A - added
        //          A$AC - added
        //           $ACG - added
        //            ACGT

        assert(get_clipping() >= gap_length);
        trim_clipping();

        sequence_ = std::string(1, '$') + sequence_;
        cigar_.data().insert(cigar_.data().begin(), Cigar::value_type{ Cigar::DELETION, 1 });
        score_ += config.gap_opening_penalty;

        if (static_cast<size_t>(gap_length) <= node_overlap) {
            // overlap is small, so add only the required dummy nods
            trim_offset();
            assert(extra_nodes >= 2);
            score_ += config.gap_opening_penalty
                + (extra_nodes - 2) * config.gap_extension_penalty;
            cigar_.data().insert(cigar_.data().begin(),
                                 Cigar::value_type{ Cigar::NODE_INSERTION, extra_nodes - 1 });
        }

        extend_query_begin(query_view_.data() - gap_length);
    }

    nodes_.insert(nodes_.begin(), extra_nodes, DeBruijnGraph::npos);

    assert(nodes_.size() == sequence_.size());
    offset_ = node_overlap;
}

// Return the string spelled by the path. This path may have disconnects (if it came)
// from a chain alignment), so this method handles that case. If there is an invalid
// edge, or if there is too long of a stretch of npos nodes, this throws a runtime error.
std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path,
                       size_t offset) {
    std::string seq;
    assert(offset < graph.get_k());

    if (path.empty())
        return seq;

    seq.reserve(path.size() + graph.get_k() - 1 - offset);

    size_t num_dummy = 0;
    if (path.front()) {
        seq += graph.get_node_sequence(path.front()).substr(offset);
    } else {
        seq += std::string(graph.get_k() - offset, '$');
        num_dummy = 1;
    }

    for (size_t i = 1; i < path.size(); ++i) {
        if (num_dummy > graph.get_k()) {
            logger->error("Too many dummy nodes\n{}", fmt::join(path, " "));
            throw std::runtime_error("");
        }

        if (path[i]) {
            if (num_dummy) {
                seq += '$';
                std::string next_seq = graph.get_node_sequence(path[i]);
                std::copy(next_seq.begin(), next_seq.end(), seq.end() - next_seq.size());
                num_dummy = 0;
            } else {
                char next = '\0';
                graph.call_outgoing_kmers(path[i - 1], [&](auto next_node, char c) {
                    if (next_node == path[i])
                        next = c;
                });

                if (!next) {
                    logger->error("Invalid edge {} {}\t{} {}",
                                  path[i - 1], path[i],
                                  graph.get_node_sequence(path[i - 1]),
                                  graph.get_node_sequence(path[i]));
                    throw std::runtime_error("");
                }

                seq += next;
            }
        } else {
            seq += '$';
            ++num_dummy;
        }
    }

    assert(seq.size() == path.size() + graph.get_k() - 1 - offset);
    return seq;
}

bool Alignment::is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config) const {
    if (empty())
        return true;

    try {
        std::string spelling = spell_path(graph, nodes_, offset_);
        if (spelling != sequence_) {
            logger->error("Stored sequence is incorrect\n{}\t{}", spelling, *this);
            return false;
        }
    } catch (const std::runtime_error&) {
        logger->error("{}", *this);
        return false;
    }

    if (!cigar_.is_valid(sequence_, query_view_)) {
        logger->error("{}", *this);
        return false;
    }

    score_t cigar_score = config ? config->score_cigar(sequence_, query_view_, cigar_) : 0;
    cigar_score += extra_penalty;
    if (config && score_ != cigar_score) {
        logger->error("Mismatch between CIGAR and score\nCigar score: {}\n{}\t{}",
                      cigar_score, query_view_, *this);
        return false;
    }

    return true;
}


AlignmentResults::AlignmentResults(std::string_view query, bool is_reverse_complement) {
    // Pad sequences for easier access in 64-bit blocks.
    // Take the max of the query size and sizeof(query_view_) to ensure that small-string
    // optimizations are disabled. Implementation taken from:
    // https://stackoverflow.com/questions/34788789/disable-stdstrings-sso
    query_.reserve(std::max(query.size(), sizeof(query_)) + 8);

    // TODO: use alphabet encoder
    // transform to upper and fix non-standard characters
    std::transform(query.begin(), query.end(), std::back_inserter(query_),
                   [](char c) { return c >= 0 ? toupper(c) : 127; });

    // fill padding with '\0'
    memset(query_.data() + query.size(), '\0', query_.capacity() - query.size());

    // set the reverse complement
    query_rc_ = query_;

    // fill padding just in case optimizations removed it
    query_rc_.reserve(query_.capacity());
    memset(query_rc_.data() + query.size(), '\0', query_rc_.capacity() - query.size());

    // reverse complement
    reverse_complement(query_rc_.begin(), query_rc_.end());
    if (is_reverse_complement)
        std::swap(query_, query_rc_);
}

} // namespace align
} // namespace graph
} // namespace mtg
