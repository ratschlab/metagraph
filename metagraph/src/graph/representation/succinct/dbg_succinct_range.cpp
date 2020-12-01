#include "dbg_succinct_range.hpp"

#include "common/seq_tools/reverse_complement.hpp"

using namespace mtg;
using namespace mtg::graph;

std::vector<boss::BOSS::TAlphabet> encode_sequence(const boss::BOSS &boss,
                                                   std::string_view sequence) {
    auto encoded = boss.encode(sequence);
    for (size_t i = 0; i < encoded.size(); ++i) {
        if (sequence[i] == boss::BOSS::kSentinel) {
            assert(encoded[i] == boss.alph_size);
            encoded[i] = boss::BOSS::kSentinelCode;
        }
    }
    return encoded;
}


DBGSuccinctRange::node_index DBGSuccinctRange
::traverse(node_index node, char next_char) const {
    if (node < offset_) {
        // normal traversal
        if (next_char != boss::BOSS::kSentinel)
            return dbg_succ_.traverse(node, next_char);

        // TODO: Use an LCS array to create a sink dummy k-mer
        return 0;
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;
    assert(offset < dbg_succ_.get_k());

    if (next_char == boss::BOSS::kSentinel) {
        if (!is_sink)
            return offset == boss_graph.get_k() ? toggle_node_sink_source(node) : 0;

        // TODO: Use an LCS array to create a sink dummy k-mer
        return 0;
    }

    if (is_sink)
        return 0;

    boss::BOSS::TAlphabet c = boss_graph.encode(next_char);
    if (c >= boss_graph.alph_size)
        return 0;

    if (!offset) {
        auto edge = boss_graph.pick_edge(last, c);
        node_index next_node = edge ? dbg_succ_.boss_to_kmer_index(edge) : 0;
        assert(!next_node
            || get_node_sequence(node).substr(1) + next_char == get_node_sequence(next_node));
        return next_node;
    }

    if (!boss_graph.tighten_range(&first, &last, c))
        return 0;

    assert(first <= last);
    assert(offset > 1 || boss_graph.succ_last(first) == last);

    std::unique_lock<std::mutex> lock(edge_pair_mutex_);
    auto it = edge_pairs_.emplace(first, last, offset - 1).first;
    node_index next_node = offset_ + ((it - edge_pairs_.begin()) * 2);
    lock.unlock();

    assert(!next_node
        || get_node_sequence(node).substr(1) + next_char == get_node_sequence(next_node));

    return next_node;
}

void DBGSuccinctRange
::call_left_tightened_ranges(node_index node,
                             boss::BOSS::edge_index first,
                             boss::BOSS::edge_index last,
                             size_t offset,
                             const std::function<void(node_index, boss::BOSS::TAlphabet)> &callback,
                             boss::BOSS::TAlphabet s) const {
    // e.g., suppose node: GTCC$$, prev_char: c
    // then traverse_back: GTCC$$ -> CGTCC$

    // In the actual implementation, the range is represented by $$GTCC, so
    // this function does $$GTCC -> $CGTCC, then sets is_sink to true

    // TODO: this can be done more efficiently with an LCS array

    std::ignore = node;
    assert(offset);
    assert(s);
#ifndef NDEBUG
    auto [edge_range, is_sink] = fetch_edge_range(node);
    assert(is_sink);
    assert(edge_range == std::tie(first, last, offset));
#endif

    const auto &boss_graph = dbg_succ_.get_boss();
    s = std::min(s, boss_graph.alph_size);

    --offset;

    auto last_seq = boss_graph.get_node_seq(last, boss_graph.get_k() - offset);
    assert(boss_graph.decode(last_seq).substr(1)
        == get_node_sequence(node - 1).substr(offset + 2));
    assert(last_seq.size() + offset == boss_graph.get_k());

    boss::BOSS::TAlphabet last_char = last_seq.front();

    if (s != boss_graph.alph_size && last_char < s)
        return;

    if (first == last) {
        if (s != boss_graph.alph_size && last_char != s)
            return;

        if (last_char == boss::BOSS::kSentinelCode)
            return;

        std::unique_lock<std::mutex> lock(edge_pair_mutex_);
        auto it = edge_pairs_.emplace(first, last, offset).first;
        size_t next_index = offset_ + (it - edge_pairs_.begin()) * 2;
        lock.unlock();

        assert(get_node_sequence(next_index).substr(offset + 1)
            == boss_graph.decode(last_seq));
        assert(fetch_edge_range(next_index + 1).second);

        callback(next_index + 1, last_char);
        return;
    }

    auto [first_char, first_minus] = boss_graph.get_minus_k_value(
        first, boss_graph.get_k() - offset + 1
    );

    if (first_char > s)
        return;

    if (first_char == last_char) {
        assert(s == boss_graph.alph_size || first_char == s);

        if (last_char == boss::BOSS::kSentinelCode)
            return;

        std::unique_lock<std::mutex> lock(edge_pair_mutex_);
        auto it = edge_pairs_.emplace(first, last, offset).first;
        node_index prev_node = offset_ + ((it - edge_pairs_.begin()) * 2) + 1;
        lock.unlock();

        assert(get_node_sequence(prev_node - 1).substr(offset + 1)
            == boss_graph.decode(last_seq));
        assert(fetch_edge_range(prev_node).second);

        callback(prev_node, last_char);
        return;
    }

    if (s == boss_graph.alph_size) {
        first_char = std::max(first_char, boss::BOSS::TAlphabet(1));
    } else {
        first_char = s;
        last_char = s;
    }

    for (boss::BOSS::TAlphabet s = first_char; s <= last_char; ++s) {
        last_seq.front() = s;
        node_index prev_node = kmer_to_node(last_seq.data(),
                                            last_seq.data() + last_seq.size(),
                                            true);

        if (prev_node) {
            assert(get_node_sequence(prev_node).rfind('$') == offset);
            assert(fetch_edge_range(prev_node + 1).second);
            callback(prev_node + 1, s);
        }
    }
}

DBGSuccinctRange::node_index DBGSuccinctRange
::traverse_back(node_index node, char prev_char) const {
    if (node < offset_) {
        // normal traversal
        if (prev_char != boss::BOSS::kSentinel)
            return dbg_succ_.traverse_back(node, prev_char);

        // TODO: Use an LCS array to create a source dummy k-mer
        return 0;
    }

    auto [edge_range, is_sink] = fetch_edge_range(node);
    if (prev_char == boss::BOSS::kSentinel) {
        if (is_sink) {
            auto [first, last, offset] = edge_range;
            return offset == dbg_succ_.get_boss().get_k()
                ? toggle_node_sink_source(node)
                : 0;
        }

        // TODO: Use an LCS array to create a source dummy k-mer
        return 0;
    }

    if (!is_sink)
        return 0;

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [first, last, offset] = edge_range;
    auto s = boss_graph.encode(prev_char);

    if (!offset) {
        return dbg_succ_.boss_to_kmer_index(
            boss_graph.pick_incoming_edge(boss_graph.bwd(last), s)
        );
    }

    node_index prev_node = 0;
    call_left_tightened_ranges(node, first, last, offset, [&](node_index prev, auto c) {
        assert(c != boss::BOSS::kSentinelCode);
        std::ignore = c;
        assert(c == s);
        assert(prev_node == 0);
        prev_node = prev;
    }, s);

    return prev_node;
}

void DBGSuccinctRange::traverse(node_index start,
                                const char *begin,
                                const char *end,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    if (start < offset_) {
        dbg_succ_.traverse(start, begin, end, callback, terminate);
    } else {
        DeBruijnGraph::traverse(start, begin, end, callback, terminate);
    }
}

void DBGSuccinctRange
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    if (terminate())
        return;

    const auto &boss_graph = dbg_succ_.get_boss();

    if (!sequence.size())
        return;

    bool is_sink = false;
    if (sequence[0] == boss::BOSS::kSentinel) {
        size_t start = sequence.rfind(boss::BOSS::kSentinel) + 1;
        sequence = std::string_view(sequence.data() + start, sequence.length() - start);
    } else if (sequence.back() == boss::BOSS::kSentinel) {
        is_sink = true;
        sequence = std::string_view(sequence.data(),
                                    sequence.find(boss::BOSS::kSentinel));
    }

    auto encoded = encode_sequence(boss_graph, sequence);

    size_t i = 0;
    size_t last_offset = dbg_succ_.get_k() + 1;
    if (sequence.size() >= dbg_succ_.get_k()) {
        // TODO: more efficient implementation
        dbg_succ_.map_to_nodes_sequentially(sequence, [&](auto node) {
            if (node) {
                assert(node
                    == dbg_succ_.kmer_to_node(std::string_view(sequence.data() + i,
                                                               dbg_succ_.get_k())));
                last_offset = 0;

            } else {
                node = kmer_to_node(encoded.data() + i,
                                    encoded.data() + i + boss_graph.get_k());

                if (node) {
                    size_t offset = get_offset(node);
                    assert(i + dbg_succ_.get_k() - offset <= sequence.size());
                    if (!last_offset || offset <= last_offset) {
                        last_offset = offset;
                    } else {
                        node = 0;
                    }

                } else {
                    last_offset = dbg_succ_.get_k() + 1;
                }
            }

            callback(node);

            ++i;
        }, terminate);
    }

    // TODO: always output the suffix matches? or only if the last k-mer was not found
    if (last_offset) {
        assert(terminate() || i + boss_graph.get_k() >= sequence.size());
        const auto *end = encoded.data() + encoded.size();

        while (!terminate() && i < sequence.size()) {
            callback(kmer_to_node(encoded.data() + i, end) + is_sink);
            ++i;
        }
    }
}

void DBGSuccinctRange::map_to_nodes(std::string_view sequence,
                                    const std::function<void(node_index)> &callback,
                                    const std::function<bool()> &terminate) const {
    if (!dbg_succ_.is_canonical_mode()) {
        map_to_nodes_sequentially(sequence, callback, terminate);
        return;
    }

    if (terminate())
        return;

    // TODO: in the reverse complement, skip those which have matched at length k
    //       in the forward mapping
    auto nodes = map_sequence_to_nodes(*this, sequence);

    std::string rev_seq(sequence.begin(), sequence.end());
    auto rev_nodes = nodes;
    reverse_complement(rev_seq, rev_nodes);
    std::reverse(rev_nodes.begin(), rev_nodes.end());

    if (nodes.size() < rev_nodes.size())
        nodes.resize(rev_nodes.size());

    if (nodes.size() > rev_nodes.size())
        rev_nodes.resize(nodes.size());

    for (size_t i = 0; i < nodes.size() && !terminate(); ++i) {
        assert(!nodes[i] == !rev_nodes[i]);
        assert(get_offset(nodes[i]) == get_offset(rev_nodes[i]));
        callback(std::min(nodes[i], rev_nodes[i]));
    }
}

size_t DBGSuccinctRange::outdegree(node_index node) const {
    if (node < offset_)
        return dbg_succ_.outdegree(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);
    if (is_sink)
        return 1;

    // TODO: better implementation: count the number of symbols in W in the range
    //       reimplement with sdsl::wt_pc::interval_symbols
    size_t count = 0;
    call_outgoing_kmers(node, [&count](auto, auto) { ++count; });
    return count;
}

bool DBGSuccinctRange::has_single_outgoing(node_index node) const {
    if (node < offset_)
        return dbg_succ_.has_single_outgoing(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;
    const auto &boss_graph = dbg_succ_.get_boss();

    if (is_sink)
        return true;

    assert(boss_graph.get_last(first));

    // check if there is only one node in the range, with only one edge
    return first == last && boss_graph.get_last(first - 1);
}

bool DBGSuccinctRange::has_multiple_outgoing(node_index node) const {
    if (node < offset_)
        return dbg_succ_.has_multiple_outgoing(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);
    if (is_sink)
        return false;

    auto [first, last, offset] = edge_range;
    assert(dbg_succ_.get_boss().get_last(first));

    // check if there are multiple nodes in the range, or if the node has
    // multiple outgoing edges
    return first != last || !dbg_succ_.get_boss().get_last(first - 1);
}

size_t DBGSuccinctRange::indegree(node_index node) const {
    if (node < offset_)
        return dbg_succ_.indegree(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);
    if (!is_sink)
        return 1;

    // TODO: more efficient
    size_t count = 0;
    call_incoming_kmers(node, [&count](auto, auto) { ++count; });
    return count;
}

bool DBGSuccinctRange::has_no_incoming(node_index node) const {
    return node < offset_ && dbg_succ_.has_no_incoming(node);
}

bool DBGSuccinctRange::has_single_incoming(node_index node) const {
    if (node < offset_)
        return dbg_succ_.has_single_incoming(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);

    if (!is_sink)
        return true;

    // TODO: more efficient
    return indegree(node) == 1;
}

DBGSuccinctRange::node_index DBGSuccinctRange
::kmer_to_node(const boss::BOSS::TAlphabet *begin,
               const boss::BOSS::TAlphabet *end,
               bool require_exact_length) const {
    const auto &boss_graph = dbg_succ_.get_boss();
    assert(begin <= end);
    assert(boss_graph.get_k() >= static_cast<size_t>(end - begin));

    if (std::find(begin, end, boss_graph.alph_size) != end)
        return 0;

    auto [first, last, seq_it] = boss_graph.index_range(begin, end);

    if (first == 0 || last == 0 || seq_it == begin
            || (require_exact_length && seq_it != end)) {
        return 0;
    } else {
        std::lock_guard<std::mutex> lock(edge_pair_mutex_);
        auto it = edge_pairs_.emplace(first, last,
                                      boss_graph.get_k() - (seq_it - begin)).first;
        return offset_ + ((it - edge_pairs_.begin()) * 2);
    }
}

DBGSuccinctRange::node_index DBGSuccinctRange::kmer_to_node(std::string_view kmer) const {
    assert(kmer.size() <= dbg_succ_.get_k());
    if (kmer.size() == dbg_succ_.get_k() && kmer.front() != boss::BOSS::kSentinel
            && kmer.back() != boss::BOSS::kSentinel) {
        return dbg_succ_.kmer_to_node(kmer);
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto encoded = encode_sequence(boss_graph, kmer);
    boss::BOSS::TAlphabet *start;
    boss::BOSS::TAlphabet *end;
    bool is_sink = false;
    if (kmer.size() < dbg_succ_.get_k()) {
        assert(kmer.front() != boss::BOSS::kSentinel);
        assert(kmer.back() != boss::BOSS::kSentinel);
        start = encoded.data();
        end = start + encoded.size();
    } else if (kmer.front() == boss::BOSS::kSentinel) {
        end = encoded.data() + encoded.size();
        start = end - (std::find(encoded.rbegin(), encoded.rend(),
                                 boss::BOSS::kSentinelCode) - encoded.rbegin());
    } else {
        assert(kmer.back() == boss::BOSS::kSentinel);
        is_sink = true;
        start = encoded.data();
        end = std::find(start, start + encoded.size(), boss::BOSS::kSentinelCode);
    }

    assert(start + boss_graph.get_k() >= end);

    node_index node = kmer_to_node(start, end, true);
    return node + (node && is_sink);
}

void DBGSuccinctRange
::call_outgoing_kmers(node_index kmer, const OutgoingEdgeCallback &callback) const {
    if (kmer < offset_) {
        dbg_succ_.call_outgoing_kmers(kmer, callback);
        return;
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [edge_range, is_sink] = fetch_edge_range(kmer);
    auto [first, last, offset] = edge_range;
    if (is_sink) {
        if (offset == boss_graph.get_k()) {
            callback(toggle_node_sink_source(kmer), boss::BOSS::kSentinel);
            return;
        }

        // TODO: Use an LCS array to create a sink dummy k-mer
        return;
    }

    if (!offset) {
        for (size_t i = boss_graph.pred_last(last - 1) + 1; i <= last; ++i) {
            auto next_node = dbg_succ_.boss_to_kmer_index(i);
            char next_char;
            if (next_node) {
                next_char = boss_graph.decode(boss_graph.get_W(i) % boss_graph.alph_size);
                if (next_char == boss::BOSS::kSentinel)
                    next_node = 0;
            }

            assert(traverse(kmer, boss_graph.decode(boss_graph.get_W(i) % boss_graph.alph_size)) == next_node);
            if (next_node)
                callback(next_node, next_char);
        }

        return;
    }

    // TODO: reimplement with sdsl::wt_pc::interval_symbols?
    for (boss::BOSS::TAlphabet c = 1; c < boss_graph.alph_size; ++c) {
        auto next_first = first;
        auto next_last = last;
        if (!boss_graph.tighten_range(&next_first, &next_last, c))
            continue;

        assert(next_first <= next_last);
        assert(offset > 1 || boss_graph.succ_last(next_first) == next_last);

        std::unique_lock<std::mutex> lock(edge_pair_mutex_);
        auto it = edge_pairs_.emplace(next_first, next_last, offset - 1).first;
        size_t next_index = offset_ + (it - edge_pairs_.begin()) * 2;
        lock.unlock();
        callback(next_index, boss_graph.decode(c));
    }
}

void DBGSuccinctRange
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    if (node < offset_) {
        dbg_succ_.adjacent_outgoing_nodes(node, callback);
        return;
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;
    if (is_sink) {
        if (offset == boss_graph.get_k()) {
            callback(toggle_node_sink_source(node));
            return;
        }

        // TODO: Use an LCS array to create a sink dummy k-mer
        return;
    }

    if (!offset) {
        for (size_t i = boss_graph.pred_last(last - 1) + 1; i <= last; ++i) {
            auto next_node = dbg_succ_.boss_to_kmer_index(i);
            assert(traverse(node, boss_graph.decode(boss_graph.get_W(i) % boss_graph.alph_size)) == next_node);
            if (next_node)
                callback(next_node);
        }

        return;
    }

    // TODO: reimplement with sdsl::wt_pc::interval_symbols?
    for (boss::BOSS::TAlphabet c = 1; c < boss_graph.alph_size; ++c) {
        auto next_first = first;
        auto next_last = last;
        if (!boss_graph.tighten_range(&next_first, &next_last, c))
            continue;

        assert(next_first <= next_last);
        assert(offset > 1 || boss_graph.succ_last(next_first) == next_last);

        std::unique_lock<std::mutex> lock(edge_pair_mutex_);
        auto it = edge_pairs_.emplace(next_first, next_last, offset - 1).first;
        size_t next_index = offset_ + (it - edge_pairs_.begin()) * 2;
        lock.unlock();
        callback(next_index);
    }
}

void DBGSuccinctRange
::call_incoming_kmers(node_index kmer, const IncomingEdgeCallback &callback) const {
    if (kmer < offset_) {
        dbg_succ_.call_incoming_kmers(kmer, callback);
        return;
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [edge_range, is_sink] = fetch_edge_range(kmer);
    auto [first, last, offset] = edge_range;
    if (!is_sink) {
        if (offset == boss_graph.get_k()) {
            callback(toggle_node_sink_source(kmer), boss::BOSS::kSentinel);
            return;
        }

        // TODO: Use an LCS array to create a source dummy k-mer
        return;
    }

    if (!offset) {
        boss_graph.call_incoming_to_target(boss_graph.bwd(last),
            boss_graph.get_node_last_value(last),
            [&](boss::BOSS::edge_index incoming_boss_edge) {
                auto prev = dbg_succ_.boss_to_kmer_index(incoming_boss_edge);
                if (prev != npos) {
                    callback(prev,
                        boss_graph.decode(
                            boss_graph.get_minus_k_value(incoming_boss_edge, boss_graph.get_k() - 1).first
                        )
                    );
                }
            }
        );

        return;
    }

    call_left_tightened_ranges(kmer, first, last, offset, [&](node_index prev, auto s) {
        assert(s != boss::BOSS::kSentinelCode);
        callback(prev, boss_graph.decode(s));
    });
}

void DBGSuccinctRange
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    if (node < offset_) {
        dbg_succ_.adjacent_incoming_nodes(node, callback);
        return;
    }

    const auto &boss_graph = dbg_succ_.get_boss();
    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;
    if (!is_sink) {
        if (offset == boss_graph.get_k()) {
            callback(toggle_node_sink_source(node));
            return;
        }

        // TODO: Use an LCS array to create a source dummy k-mer
        return;
    }

    if (!offset) {
        boss_graph.call_incoming_to_target(boss_graph.bwd(last),
            boss_graph.get_node_last_value(last),
            [&](boss::BOSS::edge_index incoming_boss_edge) {
                auto prev = dbg_succ_.boss_to_kmer_index(incoming_boss_edge);
                if (prev != npos)
                    callback(prev);
            }
        );

        return;
    }

    call_left_tightened_ranges(node, first, last, offset, [&](node_index prev, auto) {
        callback(prev);
    });
}

std::string DBGSuccinctRange::get_node_sequence(node_index node) const {
    if (node < offset_)
        return dbg_succ_.get_node_sequence(node);

    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;

    const auto &boss_graph = dbg_succ_.get_boss();

    size_t match_size = boss_graph.get_k() - offset;
    std::vector<boss::BOSS::TAlphabet> seq = boss_graph.get_node_seq(last, match_size);
    seq.resize(dbg_succ_.get_k(), boss::BOSS::kSentinelCode);
    if (!is_sink)
        std::rotate(seq.begin(), seq.begin() + match_size, seq.end());

    assert(seq.size() == dbg_succ_.get_k());

    return boss_graph.decode(seq);
}

size_t DBGSuccinctRange::get_offset(node_index node) const {
    if (node < offset_)
        return 0;

    auto [edge_range, is_sink] = fetch_edge_range(node);
    auto [first, last, offset] = edge_range;
    return offset + 1;
}

DBGSuccinctRange::EdgeDescriptor DBGSuccinctRange::fetch_edge_range(node_index node) const {
    assert(node >= offset_);
    std::lock_guard<std::mutex> lock(edge_pair_mutex_);
    assert(((node - offset_) / 2) < edge_pairs_.size());

    return std::make_pair(*(edge_pairs_.begin() + ((node - offset_) / 2)),
                          (node - offset_) % 2);
}

DBGSuccinctRange::node_index DBGSuccinctRange
::toggle_node_sink_source(node_index node) const {
    if (node < offset_)
        return node;

    size_t toggled = ((node - offset_) ^ node_index(1)) + offset_;
    assert(std::get<0>(fetch_edge_range(node)) == std::get<0>(fetch_edge_range(toggled)));
    assert(std::get<1>(fetch_edge_range(node)) != std::get<1>(fetch_edge_range(toggled)));

    return toggled;
}

bool DBGSuccinctRange::is_sink(node_index node) const {
    return node >= offset_ && ((node - offset_) & 1);
}

void DBGSuccinctRange
::reverse_complement(std::string &rev_seq, std::vector<node_index> &path) const {
#ifndef NDEBUG
    std::string orig_seq = rev_seq;
#endif

    ::reverse_complement(rev_seq.begin(), rev_seq.end());

    std::vector<node_index> rev_path = map_sequence_to_nodes(*this, rev_seq);

    if (path.size() < rev_path.size()) {
        path.resize(rev_path.size());
    } else if (path.size() > rev_path.size()) {
        rev_path.resize(path.size());
    }

    assert(rev_path.size() == path.size());

    std::vector<node_index> rev_path_init(rev_path.size());
    std::swap(rev_path, rev_path_init);

    std::vector<size_t> offsets(rev_path.size(), get_k());

    for (size_t i = 0; i < rev_path_init.size(); ++i) {
        if (!rev_path_init[i])
            continue;

        size_t offset = get_offset(rev_path_init[i]);

        // e.g.,
        // k = 6; i = 4; offset = 2
        // rev: rev_end_pos = 8
        //         --****
        // GCTAGCATCTGAGAGGGGA
        // TCCCCTCTCAGATGCTAGC
        //            ****--
        // fwd: fwd_start_pos = 11

        size_t match_length = get_k() - offset;
        size_t rev_end_pos = i + match_length;
        size_t fwd_start_pos = rev_seq.length() - rev_end_pos;
        assert(fwd_start_pos + match_length <= rev_seq.length());

        if (fwd_start_pos < rev_path.size() && offset < offsets[fwd_start_pos]) {
            offsets[fwd_start_pos] = offset;

            rev_path[fwd_start_pos] = toggle_node_sink_source(rev_path_init[i]);
#ifndef NDEBUG
            std::string tmp = get_node_sequence(rev_path[fwd_start_pos]);
            assert(tmp.front() != boss::BOSS::kSentinel);
            if (tmp.back() == boss::BOSS::kSentinel)
                tmp = tmp.substr(0, match_length);

            assert(tmp.find(boss::BOSS::kSentinel) == std::string::npos);

            ::reverse_complement(tmp.begin(), tmp.end());
            assert(tmp == std::string_view(orig_seq.data() + fwd_start_pos, match_length));
#endif
        }
    }

    std::reverse(rev_path.begin(), rev_path.end());
    std::swap(rev_path, path);
}

void DBGSuccinctRange
::call_source_nodes(const std::function<void(node_index)> &callback) const {
    const auto &boss_graph = dbg_succ_.get_boss();
    std::unique_lock<std::mutex> lock(edge_pair_mutex_);
    auto it = edge_pairs_.emplace(1, boss_graph.num_edges(), boss_graph.get_k()).first;
    node_index source = offset_ + ((it - edge_pairs_.begin()) * 2);
    lock.unlock();

    assert(get_node_sequence(source)
        == std::string(dbg_succ_.get_k(), boss::BOSS::kSentinel));

    callback(source);
}

void DBGSuccinctRange
::call_nodes_in_range(node_index node,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate) const {
    if (terminate())
        return;

    if (node < offset_) {
        callback(node);
        return;
    }

    auto edge_descriptor = fetch_edge_range(node);
    auto &[edge_range, is_sink] = edge_descriptor;
    auto &[first, last, offset] = edge_range;

    const auto &boss_graph = dbg_succ_.get_boss();

    if (is_sink) {
        assert(false);
        throw std::runtime_error("Not implemented: find all nodes with matching prefixes");
        return;
    }

    while (!terminate() && first <= last) {
        boss::BOSS::edge_index e = boss_graph.bwd(last);

        // this block is an optimization to avoid the call_incoming_to_target
        // call if only one node from the range is desired
        node_index call_val = dbg_succ_.boss_to_kmer_index(e);
        if (call_val) {
            assert(get_node_sequence(node).substr(offset + 1)
                == get_node_sequence(call_val).substr(offset + 1));
            callback(call_val);
            if (terminate())
                return;
        }

        boss_graph.call_incoming_to_target(e, boss_graph.get_node_last_value(last),
            [&](boss::BOSS::edge_index incoming_edge_idx) {
                if (!terminate() && incoming_edge_idx != e) {
                    auto call_val = dbg_succ_.boss_to_kmer_index(incoming_edge_idx);
                    if (call_val) {
#ifndef NDEBUG
                        size_t offset = std::get<2>(edge_descriptor.first);
#endif
                        assert(get_node_sequence(node).substr(offset + 1)
                            == get_node_sequence(call_val).substr(offset + 1));
                        callback(call_val);
                    }
                }
            }
        );

        if (!terminate())
            last = boss_graph.pred_last(last - 1);
    }
}
