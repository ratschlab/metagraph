#include "node_first_cache.hpp"

namespace mtg {
namespace graph {

char NodeFirstCache::get_first_char(node_index node) const {
    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

    edge_index edge = dbg_succ_->kmer_to_boss_index(node);
    assert(boss.bwd(edge) == get_parent_pair(edge).first);

    boss::BOSS::TAlphabet s = boss.get_node_last_value(get_parent_pair(edge).second);
    assert(s == boss.get_minus_k_value(edge, boss.get_k() - 1).first);

    return boss.decode(s);
}

void NodeFirstCache::call_incoming_kmers(node_index node,
                                         const IncomingEdgeCallback &callback) const {
    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

    assert(node > 0 && node <= dbg_succ_->num_nodes());

    edge_index edge = dbg_succ_->kmer_to_boss_index(node);

    edge_index bwd = 0;
    if (auto fetch = first_cache_.TryGet(edge)) {
        assert(fetch->first == boss.bwd(edge));
        bwd = fetch->first;
    } else {
        bwd = boss.bwd(edge);
    }

    boss.call_incoming_to_target(bwd, boss.get_node_last_value(edge),
        [&](boss::BOSS::edge_index incoming_boss_edge) {
            assert(boss.get_W(incoming_boss_edge) % boss.alph_size
                    == boss.get_node_last_value(edge));

            node_index prev = dbg_succ_->boss_to_kmer_index(incoming_boss_edge);
            if (prev != DeBruijnGraph::npos) {
                char c = boss.decode(boss.get_node_last_value(
                    get_parent_pair(incoming_boss_edge, edge).second
                ));
                assert(dbg_succ_->traverse_back(node, c) == prev);
                callback(prev, c);
            }
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

    if (dbg_succ != dbg_succ_) {
        if (verbose)
            std::cerr << "Graphs do not match\n";

        return false;
    }

    return true;
}

auto NodeFirstCache::get_parent_pair(edge_index edge, edge_index child_hint) const
        -> std::pair<edge_index, edge_index> {
    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

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

} // namespace graph
} // namespace mtg
