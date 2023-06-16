#include "node_rc.hpp"

#include <mutex>

#include <ips4o.hpp>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/template_utils.hpp"
#include "common/utils/file_utils.hpp"

namespace mtg {
namespace graph {

using boss::BOSS;
using edge_index = BOSS::edge_index;
using mtg::common::logger;


NodeRC::NodeRC(const DeBruijnGraph &graph) : graph_(&graph) {
    if (graph.get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only implemented for PRIMARY graphs");
        exit(1);
    }
}

void NodeRC::adjacent_outgoing_from_rc(node_index node,
                                       const std::function<void(node_index)> &callback,
                                       const NodeFirstCache *cache,
                                       const std::string &spelling_hint) const {
    //        lshift    rc
    // AGCCAT -> *AGCCA -> TGGCT*
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //-> AGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //-> TCAAGCAGAAGACGGCATACGAGATCCTCT
        const BOSS &boss = dbg_succ_->get_boss();
        edge_index rc_edge = 0;

        if (cache) {
            rc_edge = cache->get_prefix_rc(node, spelling_hint);
        } else {
            std::string rev_seq = spelling_hint.size()
                ? spelling_hint
                : dbg_succ_->get_node_sequence(node);
            rev_seq.pop_back();
            assert(rev_seq.size() == boss.get_k());

            if (rev_seq[0] == BOSS::kSentinel)
                return;

            ::reverse_complement(rev_seq.begin(), rev_seq.end());
            auto encoded = boss.encode(rev_seq);
            auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

            if (end == encoded.end()) {
                assert(edge == edge_2);
                rc_edge = edge;
            }
        }

        if (!rc_edge)
            return;

        boss.call_outgoing(rc_edge, [&](edge_index adjacent_edge) {
            assert(dbg_succ_);
            node_index prev = dbg_succ_->boss_to_kmer_index(adjacent_edge);
            if (prev != DeBruijnGraph::npos)
                callback(prev);
        });
    } else {
        // Do the checks by directly mapping the sequences of the desired k-mers.
        // For non-DBGSuccinct graphs, this should be fast enough.
        std::string rev_seq = graph_->get_node_sequence(node).substr(0, graph_->get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        rev_seq.push_back('\0');

        const auto &alphabet = graph_->alphabet();
        for (size_t c = 0; c < alphabet.size(); ++c) {
            if (alphabet[c] != boss::BOSS::kSentinel) {
                rev_seq.back() = complement(alphabet[c]);
                node_index prev = graph_->kmer_to_node(rev_seq);
                if (prev != DeBruijnGraph::npos)
                    callback(prev);
            }
        }
    }
}

void NodeRC::adjacent_incoming_from_rc(node_index node,
                                       const std::function<void(node_index)> &callback,
                                       const NodeFirstCache *cache,
                                       const std::string &spelling_hint) const {
    //        rshift    rc
    // ATGGCT -> TGGCT* -> *AGCCA
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  GAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  CTCAAGCAGAAGACGGCATACGAGATCCTC
        const BOSS &boss = dbg_succ_->get_boss();
        edge_index rc_edge = 0;

        if (cache) {
            rc_edge = cache->get_suffix_rc(node, spelling_hint);
        } else {
            std::string rev_seq = spelling_hint.size()
                ? spelling_hint.substr(1)
                : dbg_succ_->get_node_sequence(node).substr(1);

            assert(rev_seq.size() == boss.get_k());

            if (rev_seq[0] == BOSS::kSentinel)
                return;

            ::reverse_complement(rev_seq.begin(), rev_seq.end());
            auto encoded = boss.encode(rev_seq);
            auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

            if (end == encoded.end()) {
                assert(edge == edge_2);
                rc_edge = edge;
            }
        }

        if (!rc_edge)
            return;

        // rc_edge may be a dummy sink, so this won't work with graph_->call_incoming_kmers
        boss.call_incoming_to_target(
            cache ? cache->get_parent_pair(rc_edge).first : boss.bwd(rc_edge),
            boss.get_node_last_value(rc_edge),
            [&](edge_index incoming_boss_edge) {
                node_index next = dbg_succ_->boss_to_kmer_index(incoming_boss_edge);
                if (next != DeBruijnGraph::npos)
                    callback(next);
            }
        );

        return;
    }

    // Do the checks by directly mapping the sequences of the desired k-mers.
    // For non-DBGSuccinct graphs, this should be fast enough.
    std::string rev_seq = graph_->get_node_sequence(node).substr(1) + std::string(1, '\0');
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    assert(rev_seq[0] == '\0');

    const auto &alphabet = graph_->alphabet();
    for (size_t c = 0; c < alphabet.size(); ++c) {
        if (alphabet[c] != boss::BOSS::kSentinel) {
            rev_seq[0] = complement(alphabet[c]);
            node_index next = graph_->kmer_to_node(rev_seq);
            if (next != DeBruijnGraph::npos)
                callback(next);
        }
    }
}

void NodeRC::call_outgoing_from_rc(node_index node,
                                   const std::function<void(node_index, char)> &callback,
                                   const NodeFirstCache *cache,
                                   const std::string &spelling_hint) const {
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        const auto &boss = dbg_succ_->get_boss();
        adjacent_outgoing_from_rc(node, [&](node_index next) {
            char c = boss.decode(boss.get_W(dbg_succ_->kmer_to_boss_index(next)) % boss.alph_size);

            if (c == boss::BOSS::kSentinel)
                return;

            callback(next, complement(c));
        }, cache, spelling_hint);
    } else {
        adjacent_outgoing_from_rc(node, [&](node_index next) {
            char c = graph_->get_node_sequence(next).back();
            if (c == boss::BOSS::kSentinel)
                return;

            callback(next, complement(c));
        }, cache, spelling_hint);
    }
}

void NodeRC::call_incoming_from_rc(node_index node,
                                   const std::function<void(node_index, char)> &callback,
                                   const NodeFirstCache *cache,
                                   const std::string &spelling_hint) const {
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        const auto &boss = dbg_succ_->get_boss();
        adjacent_incoming_from_rc(node, [&](node_index prev) {
            char c = cache
                ? cache->get_first_char(prev)
                : boss.decode(boss.get_minus_k_value(
                    dbg_succ_->kmer_to_boss_index(prev), boss.get_k() - 1).first);

            if (c == boss::BOSS::kSentinel)
                return;

            callback(prev, complement(c));
        }, cache, spelling_hint);
    } else {
        adjacent_incoming_from_rc(node, [&](node_index prev) {
            char c = graph_->get_node_sequence(prev)[0];
            if (c == boss::BOSS::kSentinel)
                return;

            callback(prev, complement(c));
        }, cache, spelling_hint);
    }
}

bool NodeRC::is_compatible(const SequenceGraph &graph, bool) const {
    if (const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph)) {
        if (dbg == graph_)
            return true;

        logger->error("Stored graph pointer does not match");
        return false;
    }

    logger->error("Only compatible with DeBruijnGraph");
    return false;
}

} // namespace graph
} // namespace mtg
