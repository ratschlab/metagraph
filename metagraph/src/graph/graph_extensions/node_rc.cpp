#include "node_rc.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/bitmap_builder.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/serialization.hpp"

#include <filesystem>
#include <string_view>
#include <progress_bar.hpp>
#include <ips4o.hpp>
#include <mutex>

#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace graph {

static const uint64_t kBlockSize = 10'000;

NodeRC::NodeRC(const DeBruijnGraph &graph) {
    const DBGSuccinct *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);

    if (!dbg_succ) {
        common::logger->error("Only implemented for succinct graphs");
        exit(1);
    }

    if (graph.get_mode() != DeBruijnGraph::PRIMARY) {
        common::logger->error("Only implemented for PRIMARY graphs");
        exit(1);
    }

    const boss::BOSS &boss = dbg_succ->get_boss();
    static_assert(sizeof(edge_index) == sizeof(boss::BOSS::edge_index));

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
    rc_ = std::make_shared<const bit_vector_smart>(initialization_data.call_ones,
                                                   initialization_data.size,
                                                   initialization_data.num_set_bits);
    builder.reset();

    common::logger->trace("Constructing reverse complement map");
    ips4o::parallel::sort(mapping.begin(), mapping.end(), utils::LessFirst(), get_num_threads());
    std::vector<uint64_t> map;
    map.reserve(mapping.size() * 2);

    for (size_t i = 0; i < mapping.size(); ++i) {
        assert(rc_->select1(i + 1) == mapping[i].first);
        map.emplace_back(mapping[i].second.first);
        map.emplace_back(mapping[i].second.second);
    }

    mapping = decltype(mapping)();

    common::logger->trace("Compressing reverse complement map");
    mapping_ = std::make_shared<const mapping_type>(std::move(map));
}

bool NodeRC::load(const std::string &filename_base) {
    const auto rc_filename = utils::make_suffix(filename_base, kRCExtension);
    try {
        std::ifstream instream(rc_filename, std::ios::binary);
        if (!instream.good())
            return false;

        auto rc = std::make_unique<bit_vector_smart>();
        rc->load(instream);
        rc_ = std::move(rc);

        auto mapping = std::make_unique<mapping_type>();
        mapping->load(instream);
        mapping_ = std::move(mapping);
        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load graph RC from file "
                  << rc_filename << std::endl;
        return false;
    }
}

void NodeRC::serialize(const std::string &filename_base) const {
    const auto fname = utils::make_suffix(filename_base, kRCExtension);

    std::ofstream outstream(fname, std::ios::binary);
    rc_->serialize(outstream);
    mapping_->serialize(outstream);
}

bool NodeRC::is_compatible(const SequenceGraph &graph, bool verbose) const {
    if (graph.max_index() + 1 != rc_->size()) {
        if (verbose)
            std::cerr << "ERROR: RC file does not match number of nodes in graph"
                      << std::endl;
        return false;
    }

    if (rc_->num_set_bits() * 2 == mapping_->size()) {
        if (verbose)
            std::cerr << "ERROR: RC file contains the wrong mapping"
                      << std::endl;
    }

    return true;
}

} // namespace graph
} // namespace mtg
