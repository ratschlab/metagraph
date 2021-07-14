#include "aligner_alignment.hpp"

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/rc_dbg.hpp"


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
void Alignment<NodeType>::append(Alignment&& other) {
    assert(query_.data() + query_.size() + other.get_clipping() == other.query_.data());
    assert(orientation_ == other.orientation_);

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;

    cigar_.append(std::move(other.cigar_));
    query_ = { query_.data(), other.query_.size() + (other.query_.data() - query_.data()) };
}

template <typename NodeType>
size_t Alignment<NodeType>::trim_offset() {
    assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos) == nodes_.end()
            && "chains not supported");

    if (!offset_ || nodes_.size() <= 1)
        return 0;

    size_t trim = std::min(offset_, nodes_.size() - 1);
    offset_ -= trim;
    nodes_.erase(nodes_.begin(), nodes_.begin() + trim);
    return trim;
}

template <typename NodeType>
size_t Alignment<NodeType>
::trim_query_prefix(size_t n, const DeBruijnGraph &graph, const DBGAlignerConfig &config) {
    assert(!offset_);

    size_t clipping = get_clipping() + n;

    auto it = cigar_.begin() + static_cast<bool>(clipping);
    size_t cigar_offset = 0;

    auto s_it = sequence_.begin();
    auto node_it = nodes_.begin();

    size_t offset_cutoff = graph.get_k() - 1;

    auto consume_ref = [&]() {
        assert(s_it != sequence_.end());
        ++s_it;
        if (offset_ < offset_cutoff) {
            ++offset_;
        } else if (node_it + 1 < nodes_.end()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n) {
        if (it == cigar_.end()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.end());
                score_ -= config.get_row(query_[0])[*s_it];
                query_.remove_prefix(1);
                --n;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_.remove_prefix(1);
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

    nodes_.erase(nodes_.begin(), node_it);
    sequence_.erase(sequence_.begin(), s_it);
    it->second -= cigar_offset;
    cigar_.erase(cigar_.begin(), it);

    assert(is_valid(graph, &config));

    return cigar_offset;
}

template <typename NodeType>
void Alignment<NodeType>::reverse_complement(const DeBruijnGraph &graph,
                                             std::string_view query_rev_comp) {
    assert(query_.size() + get_end_clipping() == query_rev_comp.size() - get_clipping());

    trim_offset();
    assert(!offset_ || nodes_.size() == 1);

    if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(&graph)) {
        if (offset_) {
            *this = Alignment<NodeType>();
        } else {
            std::reverse(cigar_.begin(), cigar_.end());
            std::reverse(nodes_.begin(), nodes_.end());
            ::reverse_complement(sequence_.begin(), sequence_.end());
            assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

            orientation_ = !orientation_;
            query_ = { query_rev_comp.data() + get_clipping(),
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

            const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
            const auto &dbg_succ = dynamic_cast<const DBGSuccinct&>(graph.get_base_graph());

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
                NodeType next_node = 0;
                char last_char;
                canonical->call_outgoing_kmers(nodes_[0], [&](NodeType next, char c) {
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
    if (sequence_.find("$") != std::string::npos
            || std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos) != nodes_.end()) {
        throw std::runtime_error("JSON output for chains not supported");
    }

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
            graph.call_outgoing_kmers(*(nodes_.rbegin() + 1),
                                      [&](auto node, char c) {
                if (node == nodes_.back())
                    sequence_ += c;
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
void Alignment<NodeType>::insert_gap_prefix(ssize_t gap_length,
                                            const DeBruijnGraph &graph,
                                            const DBGAlignerConfig &config) {
    size_t extra_nodes = graph.get_k();

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
            cigar_.insert(cigar_.begin(),
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
        cigar_.insert(cigar_.begin(), Cigar::value_type{ Cigar::DELETION, 1 });
        score_ += config.gap_opening_penalty;

        if (static_cast<size_t>(gap_length) < graph.get_k()) {
            // overlap is small, so add only the required dummy nods
            trim_offset();
            assert(extra_nodes >= 2);
            score_ += config.gap_opening_penalty
                + (extra_nodes - 2) * config.gap_extension_penalty;
            cigar_.insert(cigar_.begin(),
                          Cigar::value_type{ Cigar::NODE_INSERTION, extra_nodes - 1 });
        }

        extend_query_begin(query_.data() - gap_length);
    }

    nodes_.insert(nodes_.begin(), extra_nodes, DeBruijnGraph::npos);

    assert(nodes_.size() == sequence_.size());
    offset_ = graph.get_k() - 1;
}


bool spell_path(const DeBruijnGraph &graph,
                const std::vector<DeBruijnGraph::node_index> &path,
                std::string &seq,
                size_t offset) {
    assert(offset < graph.get_k());

    if (path.empty())
        return "";

    seq.clear();
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
            std::string out_msg;
            for (DeBruijnGraph::node_index node : path) {
                out_msg += fmt::format("{} ", node);
            }
            logger->error("Too many dummy nodes\n{}", out_msg);
            return false;
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
                    return false;
                }

                seq += next;
            }
        } else {
            seq += '$';
            ++num_dummy;
        }
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
template class QueryAlignment<>;

} // namespace align
} // namespace graph
} // namespace mtg
