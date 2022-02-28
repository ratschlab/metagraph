#include "node_rc.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/bitmap_builder.hpp"
#include "common/utils/template_utils.hpp"

#include <progress_bar.hpp>
#include <ips4o.hpp>
#include <mutex>

namespace mtg {
namespace graph {

using boss::BOSS;
using edge_index = BOSS::edge_index;

template <class Indicator, class Mapping>
NodeRC<Indicator, Mapping>::NodeRC(const DBGSuccinct &graph)
      : graph_(&graph) {
    if (graph.get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only implemented for PRIMARY graphs");
        exit(1);
    }

    const BOSS &boss = graph.get_boss();

    std::unique_ptr<bitmap_builder> builder;
    builder = std::make_unique<bitmap_builder_set>(graph.max_index() + 1, 0,
                                                   graph.num_nodes() * 0.06);

    std::mutex mu;

    std::vector<std::pair<node_index, std::pair<edge_index, edge_index>>> mapping;
    graph.call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
        std::string rev_seq = seq;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss.encode(rev_seq);

        // initialize
        auto it = encoded.end() - boss.get_k();
        auto [edge, edge_2, end] = boss.index_range(it, it + boss.get_k());
        assert(end != it + boss.get_k() || edge == edge_2);

        for (size_t i = 0; i < path.size(); ++i) {
            // prefix
            size_t map_i = 0;
            bool found = false;
            if (end == it + boss.get_k()) {
                std::lock_guard<std::mutex> lock(mu);
                assert(edge);
                map_i = mapping.size();
                mapping.emplace_back(path[i], std::make_pair(edge, 0));
                found = true;
                builder->add_one(path[i]);
            }

            // suffix
            --it;
            std::tie(edge, edge_2, end) = boss.index_range(it, it + boss.get_k());
            assert(end != it + boss.get_k() || edge == edge_2);

            if (end == it + boss.get_k()) {
                std::lock_guard<std::mutex> lock(mu);
                assert(edge);
                if (found) {
                    assert(mapping[map_i].first == path[i]);
                    assert(mapping[map_i].second.first);
                    assert(!mapping[map_i].second.second);
                    mapping[map_i].second.second = edge;
                } else {
                    builder->add_one(path[i]);
                    mapping.emplace_back(path[i], std::make_pair(0, edge));
                }
            }
        }
        assert(it == encoded.begin());
    }, get_num_threads());

    common::logger->trace("Found {} / {} ({:.2f}%) nodes with reverse complements",
                          mapping.size(), graph.num_nodes(),
                          100.0 * mapping.size() / graph.num_nodes());

    common::logger->trace("Generating reverse complement indicator vector");
    auto initialization_data = builder->get_initialization_data();
    assert(initialization_data.num_set_bits == mapping.size());
    rc_ = Indicator(initialization_data.call_ones,
                    initialization_data.size,
                    initialization_data.num_set_bits);
    builder.reset();

    common::logger->trace("Constructing reverse complement map");
    ips4o::parallel::sort(mapping.begin(), mapping.end(), utils::LessFirst(),
                          get_num_threads());
    std::vector<uint64_t> map;
    map.reserve(mapping.size() * 2);

    for (size_t i = 0; i < mapping.size(); ++i) {
        assert(rc_.select1(i + 1) == mapping[i].first);
        map.emplace_back(mapping[i].second.first);
        map.emplace_back(mapping[i].second.second);
    }

    mapping = decltype(mapping)();

    common::logger->trace("Compressing reverse complement map");
    mapping_ = { std::move(map) };
}

template <class Indicator, class Mapping>
void NodeRC<Indicator, Mapping>
::call_outgoing_nodes_from_rc(node_index node,
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

template <class Indicator, class Mapping>
void NodeRC<Indicator, Mapping>
::call_incoming_nodes_from_rc(node_index node,
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

template <class Indicator, class Mapping>
bool NodeRC<Indicator, Mapping>::load(const std::string &filename_base) {
    const auto rc_filename = utils::make_suffix(filename_base, kRCExtension);
    try {
        std::ifstream instream(rc_filename, std::ios::binary);
        if (!instream.good())
            return false;

        rc_.load(instream);
        mapping_.load(instream);
        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load graph RC from file " << rc_filename << std::endl;
        return false;
    }
}

template <class Indicator, class Mapping>
void NodeRC<Indicator, Mapping>::serialize(const std::string &filename_base) const {
    if (!rc_.size())
        std::cerr << "WARNING: NodeRC was initialized with set_graph, so nothing to serialize." << std::endl;

    const auto fname = utils::make_suffix(filename_base, kRCExtension);

    std::ofstream outstream(fname, std::ios::binary);
    rc_.serialize(outstream);
    mapping_.serialize(outstream);
}

template <class Indicator, class Mapping>
bool NodeRC<Indicator, Mapping>::is_compatible(const SequenceGraph &graph,
                                               bool verbose) const {
    if (!rc_.size()) {
        if (const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph))
            return *dbg == *graph_;

        if (verbose)
            std::cerr << "ERROR: only compatible with DeBruijnGraph" << std::endl;

        return false;
    }

    if (graph.max_index() + 1 != rc_.size()) {
        if (verbose) {
            std::cerr << "ERROR: RC file does not match number of nodes in graph"
                      << std::endl;
        }
        return false;
    }

    if (rc_.num_set_bits() * 2 != mapping_.size()) {
        if (verbose)
            std::cerr << "ERROR: RC file contains the wrong mapping" << std::endl;

        return false;
    }

    return true;
}

template class NodeRC<>;

} // namespace graph
} // namespace mtg
