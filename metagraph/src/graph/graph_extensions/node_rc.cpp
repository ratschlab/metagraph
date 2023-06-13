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


NodeRC::NodeRC(const DeBruijnGraph &graph, bool construct_index) : graph_(&graph) {
    if (graph.get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only implemented for PRIMARY graphs");
        exit(1);
    }

    if (!construct_index)
        return;

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!dbg_succ)
        return;

    const BOSS &boss = dbg_succ->get_boss();

    std::mutex mu;

    std::vector<node_index> rc_nodes_prefix;
    std::vector<node_index> rc_nodes_suffix;
    graph.call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
        assert(seq.size() == path.size() + boss.get_k());
        std::string rev_seq = seq;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());

        // map each (k-1)-mer from the reverse complement to BOSS nodes
        auto encoded = boss.encode(rev_seq);
        std::vector<edge_index> rc_edges;
        rc_edges.reserve(path.size() + 1);
        for (auto it = encoded.rbegin() + boss.get_k(); it <= encoded.rend(); ++it) {
            auto [edge, edge_2, end] = boss.index_range(it.base(), it.base() + boss.get_k());
            assert(end != it.base() + boss.get_k() || edge == edge_2);
            rc_edges.emplace_back(end == it.base() + boss.get_k() ? edge : 0);
        }
        assert(rc_edges.size() == path.size() + 1);

        // For each contig node, determine if the prefix or the suffix of its reverse
        // complement maps to a BOSS node. If so, then record that node
        auto it = rc_edges.begin();
        for (size_t i = 0; i < path.size(); ++i, ++it) {
            if (*it) {
                std::lock_guard<std::mutex> lock(mu);
                rc_nodes_prefix.emplace_back(path[i]);
            }

            if (*(it + 1)) {
                std::lock_guard<std::mutex> lock(mu);
                rc_nodes_suffix.emplace_back(path[i]);
            }
        }
    }, get_num_threads());

    logger->trace("Found {} / {} ({:.2f}%) nodes with reverse complements",
                  rc_nodes_prefix.size() + rc_nodes_suffix.size(), graph.num_nodes(),
                  100.0 * (rc_nodes_prefix.size() + rc_nodes_suffix.size()) / graph.num_nodes());

    // sort by node ID, then construct the binary indicator vector
    logger->trace("Constructing reverse complement map");
    ips4o::parallel::sort(rc_nodes_prefix.begin(), rc_nodes_prefix.end(),
                          std::less<node_index>(), get_num_threads());
    ips4o::parallel::sort(rc_nodes_suffix.begin(), rc_nodes_suffix.end(),
                          std::less<node_index>(), get_num_threads());

    rc_prefix_ = Indicator([&](const auto &callback) {
                               std::for_each(rc_nodes_prefix.begin(),
                                             rc_nodes_prefix.end(),
                                             callback);
                           },
                           graph.max_index() + 1,
                           rc_nodes_prefix.size());
    rc_suffix_ = Indicator([&](const auto &callback) {
                               std::for_each(rc_nodes_suffix.begin(),
                                             rc_nodes_suffix.end(),
                                             callback);
                           },
                           graph.max_index() + 1,
                           rc_nodes_suffix.size());
}

void NodeRC::adjacent_outgoing_from_rc(node_index node,
                                       const std::function<void(node_index)> &callback) const {
    //        lshift    rc
    // AGCCAT -> *AGCCA -> TGGCT*
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //-> AGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //-> TCAAGCAGAAGACGGCATACGAGATCCTCT
        const BOSS &boss = dbg_succ_->get_boss();
        edge_index rc_edge = 0;

        if (rc_prefix_.size() && !rc_prefix_[node])
            return;

        std::string rev_seq;
        if (auto first_cache = dbg_succ_->get_extension<NodeFirstCache>()) {
            rev_seq = first_cache->get_node_sequence(node);
        } else {
            rev_seq = dbg_succ_->get_node_sequence(node);
        }
        rev_seq.pop_back();
        assert(rev_seq.size() == boss.get_k());

        assert(!rc_prefix_.size() || rev_seq[0] != BOSS::kSentinel);
        if (rev_seq[0] == BOSS::kSentinel)
            return;

        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss.encode(rev_seq);
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

        assert(!rc_prefix_.size() || end == encoded.end());
        if (end == encoded.end()) {
            assert(edge == edge_2);
            rc_edge = edge;
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
                                       const std::function<void(node_index)> &callback) const {
    //        rshift    rc
    // ATGGCT -> TGGCT* -> *AGCCA
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  GAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  CTCAAGCAGAAGACGGCATACGAGATCCTC
        const BOSS &boss = dbg_succ_->get_boss();
        edge_index rc_edge = 0;

        if (rc_suffix_.size() && !rc_suffix_[node])
            return;

        std::string rev_seq;
        auto first_cache = dbg_succ_->get_extension<NodeFirstCache>();
        if (first_cache) {
            rev_seq = first_cache->get_node_sequence(node).substr(1);
        } else {
            rev_seq = dbg_succ_->get_node_sequence(node).substr(1);
        }
        assert(rev_seq.size() == boss.get_k());

        assert(!rc_suffix_.size() || rev_seq[0] != BOSS::kSentinel);
        if (rev_seq[0] == BOSS::kSentinel)
            return;

        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss.encode(rev_seq);
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());

        assert(!rc_suffix_.size() || end == encoded.end());
        if (end == encoded.end()) {
            assert(edge == edge_2);
            rc_edge = edge;
        }

        if (!rc_edge)
            return;

        // rc_edge may be a dummy sink, so this won't work with graph_->call_incoming_kmers
        boss.call_incoming_to_target(
            first_cache ? first_cache->get_parent_pair(rc_edge).first : boss.bwd(rc_edge),
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
                                   const std::function<void(node_index, char)> &callback) const {
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        const auto &boss = dbg_succ_->get_boss();
        adjacent_outgoing_from_rc(node, [&](node_index next) {
            char c = boss.decode(boss.get_W(dbg_succ_->kmer_to_boss_index(next)) % boss.alph_size);

            if (c == boss::BOSS::kSentinel)
                return;

            callback(next, complement(c));
        });
    } else {
        adjacent_outgoing_from_rc(node, [&](node_index next) {
            char c = graph_->get_node_sequence(next).back();
            if (c == boss::BOSS::kSentinel)
                return;

            callback(next, complement(c));
        });
    }
}

void NodeRC::call_incoming_from_rc(node_index node,
                                   const std::function<void(node_index, char)> &callback) const {
    assert(graph_);

    if (const auto *dbg_succ_ = dynamic_cast<const DBGSuccinct*>(graph_)) {
        const auto &boss = dbg_succ_->get_boss();
        auto first_cache = dbg_succ_->get_extension<NodeFirstCache>();
        adjacent_incoming_from_rc(node, [&](node_index prev) {
            char c = first_cache
                ? first_cache->get_first_char(prev)
                : boss.decode(boss.get_minus_k_value(
                    dbg_succ_->kmer_to_boss_index(prev), boss.get_k() - 1).first);

            if (c == boss::BOSS::kSentinel)
                return;

            callback(prev, complement(c));
        });
    } else {
        adjacent_incoming_from_rc(node, [&](node_index prev) {
            char c = graph_->get_node_sequence(prev)[0];
            if (c == boss::BOSS::kSentinel)
                return;

            callback(prev, complement(c));
        });
    }
}

bool NodeRC::load(const std::string &filename_base) {
    const auto rc_filename = utils::make_suffix(filename_base, kRCExtension);
    try {
        std::unique_ptr<std::ifstream> in = utils::open_ifstream(rc_filename);
        if (!in->good())
            return false;

        rc_prefix_.load(*in);
        rc_suffix_.load(*in);
        return true;

    } catch (...) {
        return false;
    }
}

void NodeRC::serialize(const std::string &filename_base) const {
    if (!rc_prefix_.size() || !rc_suffix_.size())
        logger->warn("NodeRC was initialized with set_graph, so nothing to serialize.");

    const auto fname = utils::make_suffix(filename_base, kRCExtension);

    std::ofstream out = utils::open_new_ofstream(fname);
    rc_prefix_.serialize(out);
    rc_suffix_.serialize(out);
}

bool NodeRC::is_compatible(const SequenceGraph &graph, bool) const {
    if (!rc_prefix_.size() || !rc_suffix_.size()) {
        if (const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph)) {
            if (dbg == graph_)
                return true;

            logger->error("Stored graph pointer does not match");
            return false;
        }

        logger->error("Only compatible with DeBruijnGraph");
        return false;
    }

    if (graph.max_index() + 1 != rc_prefix_.size()) {
        logger->error("RC file does not match number of nodes in graph");
        return false;
    }

    if (graph.max_index() + 1 != rc_suffix_.size()) {
        logger->error("RC file does not match number of nodes in graph");
        return false;
    }

    return true;
}

} // namespace graph
} // namespace mtg
