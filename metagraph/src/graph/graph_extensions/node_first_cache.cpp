#include "node_first_cache.hpp"

#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {

char NodeFirstCache::get_first_char(node_index node, node_index child_hint) const {
    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

    edge_index edge = dbg_succ_->kmer_to_boss_index(node);
    assert(boss.bwd(edge) == get_parent_pair(edge).first);

    boss::BOSS::TAlphabet s = boss.get_node_last_value(get_parent_pair(edge,
        child_hint ? dbg_succ_->kmer_to_boss_index(child_hint)
                   : boss.fwd(edge)
    ).second);
    assert(s == boss.get_minus_k_value(edge, boss.get_k() - 1).first);

    return boss.decode(s);
}

std::string NodeFirstCache::get_node_sequence(node_index node) const {
    assert(dbg_succ_);
    assert(node > 0 && node <= dbg_succ_->num_nodes());

    boss::BOSS::edge_index x = dbg_succ_->kmer_to_boss_index(node);

    if (!first_cache_.TryGet(x))
        return dbg_succ_->get_node_sequence(node);

    const auto &boss = dbg_succ_->get_boss();

    std::vector<boss::BOSS::TAlphabet> ret(boss.get_k() + 1);
    ret.back() = boss.get_W(x) % boss.alph_size;

    size_t i = boss.get_k();
    ret[--i] = boss.get_node_last_value(x);

    boss::BOSS::edge_index last_edge = 0;
    while (i > 0) {
        boss::BOSS::edge_index bwd = get_parent_pair(x, last_edge).first;
        last_edge = x;
        x = bwd;
        ret[--i] = boss.get_node_last_value(x);
    }

    return boss.decode(ret);
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
            if (prev != DeBruijnGraph::npos)
                callback(prev, get_first_char(prev, node));
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

auto NodeFirstCache
::get_prefix_rc(node_index node, const std::string &spelling_hint) const -> edge_index {
    if (auto fetch = prefix_rc_cache_.TryGet(node))
        return *fetch;

    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

    std::string rev_seq = spelling_hint.size()
        ? spelling_hint
        : get_node_sequence(node);
    rev_seq.pop_back();
    assert(rev_seq.size() == boss.get_k());

    if (rev_seq[0] == boss::BOSS::kSentinel) {
        prefix_rc_cache_.Put(node, 0);
        return 0;
    }

    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    auto encoded = boss.encode(rev_seq);
    auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

    if (end != encoded.end()) {
        prefix_rc_cache_.Put(node, 0);
        return 0;
    }

    prefix_rc_cache_.Put(node, edge);
    return edge;
}

auto NodeFirstCache
::get_suffix_rc(node_index node, const std::string &spelling_hint) const -> edge_index {
    if (auto fetch = suffix_rc_cache_.TryGet(node))
        return *fetch;

    assert(dbg_succ_);
    const boss::BOSS &boss = dbg_succ_->get_boss();

    std::string rev_seq = spelling_hint.size()
        ? spelling_hint.substr(1)
        : get_node_sequence(node).substr(1);
    assert(rev_seq.size() == boss.get_k());

    if (rev_seq[0] == boss::BOSS::kSentinel) {
        suffix_rc_cache_.Put(node, 0);
        return 0;
    }

    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    auto encoded = boss.encode(rev_seq);
    auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

    if (end != encoded.end()) {
        suffix_rc_cache_.Put(node, 0);
        return 0;
    }

    suffix_rc_cache_.Put(node, edge);
    return edge;
}

} // namespace graph
} // namespace mtg
