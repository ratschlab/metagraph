#include "node_first_cache.hpp"

#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {

using namespace boss;

char NodeFirstCache::get_first_char(edge_index edge, edge_index child_hint) const {
    const BOSS &boss = dbg_succ_.get_boss();

    if (!cache_size_)
        return boss.decode(boss.get_minus_k_value(edge, boss.get_k() - 1).first);

    assert(boss.bwd(edge) == get_parent_pair(edge).first);

    BOSS::TAlphabet s = boss.get_node_last_value(get_parent_pair(edge, child_hint).second);
    assert(s == boss.get_minus_k_value(edge, boss.get_k() - 1).first);

    return boss.decode(s);
}

void NodeFirstCache::call_incoming_edges(edge_index edge,
                                         const BOSS::Call<edge_index> &callback) const {
    const BOSS &boss = dbg_succ_.get_boss();
    edge_index bwd = cache_size_ ? get_parent_pair(edge).first : boss.bwd(edge);

    boss.call_incoming_to_target(bwd, boss.get_node_last_value(edge),
        [&](BOSS::edge_index prev) {
            assert(boss.get_W(prev) % boss.alph_size == boss.get_node_last_value(edge));
            callback(prev);
        }
    );
}

void NodeFirstCache::call_incoming_kmers(node_index node,
                                         const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= dbg_succ_.num_nodes());

    edge_index edge = dbg_succ_.kmer_to_boss_index(node);

    call_incoming_edges(edge,
        [&](edge_index prev_edge) {
            node_index prev = dbg_succ_.boss_to_kmer_index(prev_edge);
            if (prev != DeBruijnGraph::npos)
                callback(prev, get_first_char(prev_edge, edge));
        }
    );
}

bool NodeFirstCache::is_compatible(const SequenceGraph &graph, bool verbose) const {
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!dbg_succ) {
        if (verbose)
            std::cerr << "Graph is not DBGSuccinct\n";

        return false;
    }

    if (dbg_succ != &dbg_succ_) {
        if (verbose)
            std::cerr << "Graphs do not match\n";

        return false;
    }

    return true;
}

auto NodeFirstCache::get_parent_pair(edge_index edge, edge_index child_hint) const
        -> std::pair<edge_index, edge_index> {
    assert(cache_size_);

    const BOSS &boss = dbg_succ_.get_boss();

    if (boss.get_k() == 1)
        return std::make_pair(boss.bwd(edge), edge);

    if (auto fetch = first_cache_.TryGet(edge)) {
        assert(fetch->first == boss.bwd(edge));
        return *fetch;
    }

    // if a child of edge has been cached already, check if we can compute
    // get_parent_pair with only two bwd calls
    if (child_hint) {
        if (auto fetch = first_cache_.TryGet(child_hint)) {
            const auto &[cur, first] = *fetch;
            if (cur == edge) {
                auto ret_val = std::make_pair(boss.bwd(edge), boss.bwd(first));
                assert(boss.get_node_last_value(ret_val.second)
                    == boss.get_minus_k_value(edge, boss.get_k() - 1).first);
                first_cache_.Put(edge, ret_val);
                return ret_val;
            }
        }
    }

    // we have to take all k - 1 steps
    assert(boss.get_k() >= 2);
    edge_index i = boss.bwd(edge);
    edge_index parent = i;
    for (size_t k = boss.get_k() - 2; k > 0; --k) {
        i = boss.bwd(i);
    }

    assert(std::make_pair(boss.get_node_last_value(i), boss.bwd(i))
        == boss.get_minus_k_value(edge, boss.get_k() - 1));
    auto ret_val = std::make_pair(parent, i);
    first_cache_.Put(edge, ret_val);

    return ret_val;
}

auto NodeFirstCache
::get_prefix_rc(edge_index edge, const std::string &spelling) const -> edge_index {
    if (auto fetch = prefix_rc_cache_.TryGet(edge))
        return *fetch;

    const BOSS &boss = dbg_succ_.get_boss();

    assert(spelling.size() == dbg_succ_.get_k());
    std::string rev_seq = spelling;
    rev_seq.pop_back();

    if (rev_seq[0] == BOSS::kSentinel) {
        prefix_rc_cache_.Put(edge, 0);
        return 0;
    }

    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    auto encoded = boss.encode(rev_seq);
    auto [rc_edge_1, rc_edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

    if (end != encoded.end()) {
        prefix_rc_cache_.Put(edge, 0);
        return 0;
    }

    prefix_rc_cache_.Put(edge, rc_edge_2);
    return rc_edge_2;
}

auto NodeFirstCache
::get_suffix_rc(edge_index edge, const std::string &spelling) const -> edge_index {
    if (auto fetch = suffix_rc_cache_.TryGet(edge))
        return *fetch;

    const BOSS &boss = dbg_succ_.get_boss();

    assert(spelling.size() == dbg_succ_.get_k());
    std::string rev_seq = spelling.substr(1);

    if (rev_seq[0] == BOSS::kSentinel) {
        suffix_rc_cache_.Put(edge, 0);
        return 0;
    }

    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    auto encoded = boss.encode(rev_seq);
    auto [rc_edge_1, rc_edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

    if (end != encoded.end()) {
        suffix_rc_cache_.Put(edge, 0);
        return 0;
    }

    suffix_rc_cache_.Put(edge, rc_edge_1);
    return rc_edge_1;
}

} // namespace graph
} // namespace mtg
