#include "node_rc.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/template_utils.hpp"

#include <ips4o.hpp>
#include <mutex>

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

    std::vector<std::pair<node_index, std::pair<edge_index, edge_index>>> mapping;
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
        // complement maps to a BOSS node. If so, then record the mapping
        auto it = rc_edges.begin();
        for (size_t i = 0; i < path.size(); ++i, ++it) {
            if (*it || *(it + 1)) {
                std::lock_guard<std::mutex> lock(mu);
                mapping.emplace_back(path[i], std::make_pair(*it, *(it + 1)));
            }
        }
    }, get_num_threads());

    logger->trace("Found {} / {} ({:.2f}%) nodes with reverse complements",
                  mapping.size(), graph.num_nodes(),
                  100.0 * mapping.size() / graph.num_nodes());

    // sort by node ID, then construct the binary indicator vector at the compressed
    // node-to-rc mapping
    logger->trace("Constructing reverse complement map");
    ips4o::parallel::sort(mapping.begin(), mapping.end(), utils::LessFirst(),
                          get_num_threads());

    std::vector<uint64_t> map;
    map.reserve(mapping.size() * 2);

    rc_ = Indicator([&](const auto &callback) {
                        for (const auto &[node, rc] : mapping) {
                            callback(node);
                            map.emplace_back(rc.first);
                            map.emplace_back(rc.second);
                        }
                    },
                    graph.max_index() + 1,
                    mapping.size());

    mapping = decltype(mapping)();

    logger->trace("Compressing reverse complement map");
    mapping_ = { std::move(map) };
}

void NodeRC::call_outgoing_from_rc(node_index node,
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

        if (rc_.size()) {
            if (auto rank = rc_.conditional_rank1(node))
                rc_edge = mapping_[(rank - 1) * 2];
        } else {
            std::string rev_seq = dbg_succ_->get_node_sequence(node).substr(0, boss.get_k());
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

        return;
    }

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

void NodeRC::call_incoming_from_rc(node_index node,
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

        if (rc_.size()) {
            if (auto rank = rc_.conditional_rank1(node))
                rc_edge = mapping_[(rank - 1) * 2 + 1];
        } else {
            std::string rev_seq = dbg_succ_->get_node_sequence(node).substr(1);
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
        boss.call_incoming_to_target(boss.bwd(rc_edge), boss.get_node_last_value(rc_edge),
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

bool NodeRC::load(const std::string &filename_base) {
    const auto rc_filename = utils::make_suffix(filename_base, kRCExtension);
    try {
        std::ifstream instream(rc_filename, std::ios::binary);
        if (!instream.good())
            return false;

        rc_.load(instream);
        mapping_.load(instream);
        return true;

    } catch (...) {
        logger->error("Cannot load graph RC from file {}", rc_filename);
        return false;
    }
}

void NodeRC::serialize(const std::string &filename_base) const {
    if (!rc_.size())
        logger->warn("NodeRC was initialized with set_graph, so nothing to serialize.");

    const auto fname = utils::make_suffix(filename_base, kRCExtension);

    std::ofstream outstream(fname, std::ios::binary);
    rc_.serialize(outstream);
    mapping_.serialize(outstream);
}

bool NodeRC::is_compatible(const SequenceGraph &graph, bool) const {
    if (!rc_.size()) {
        if (const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph)) {
            if (dbg == graph_)
                return true;

            logger->error("Stored graph pointer does not match");
            return false;
        }

        logger->error("Only compatible with DeBruijnGraph");
        return false;
    }

    if (graph.max_index() + 1 != rc_.size()) {
        logger->error("RC file does not match number of nodes in graph");
        return false;
    }

    if (rc_.num_set_bits() * 2 != mapping_.size()) {
        logger->error("RC file contains the wrong mapping");
        return false;
    }

    return true;
}

} // namespace graph
} // namespace mtg
