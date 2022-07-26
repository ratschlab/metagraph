#ifndef __UNITIGS_HPP__
#define __UNITIGS_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/aligner_seeder_methods.hpp"


namespace mtg {
namespace graph {
namespace align {

class Unitigs : public SequenceGraph::GraphExtension {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef std::function<void(std::function<void(const std::vector<node_index>&)>)> PathGenerator;

    Unitigs(const DeBruijnGraph &graph) : graph_(graph) {}

    Unitigs(const DeBruijnGraph &graph, const PathGenerator &generate_paths)
          : Unitigs(graph) {
        const DeBruijnGraph *base_dbg = &graph_;
        if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_dbg))
            base_dbg = &canonical->get_graph();

        unitigs_ = sdsl::int_vector<>(base_dbg->num_nodes() + 1, 0);

        size_t i = 1;
        generate_paths([&](const std::vector<node_index> &path) {
            for (node_index node : path) {
                unitigs_[node] = i;
            }
            ++i;
        });

        sdsl::util::bit_compress(unitigs_);
    }

    bool load(const std::string &filename_base) {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ifstream fin(fname, std::ios::binary);
        if (!fin.good())
            return false;

        try {
            unitigs_.load(fin);
            return true;
        } catch (...) {
            return false;
        }
    }

    void serialize(const std::string &filename_base) const {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ofstream fout(fname, std::ios::binary);
        unitigs_.serialize(fout);
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    std::vector<Seed> cluster_and_filter_seeds(const std::vector<Seed> &seeds,
                                               size_t min_seed_length) const {
        const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph_);

        common::logger->trace("Clustering: old seed count: {}", seeds.size());
        tsl::hopscotch_map<size_t, std::vector<size_t>> unitig_to_bucket;
        for (size_t i = 0; i < seeds.size(); ++i) {
            node_index node = seeds[i].get_nodes().back();
            if (canonical)
                node = canonical->get_base_node(node);

            size_t unitig_id = unitigs_[node];
            unitig_to_bucket[unitig_id].emplace_back(i);
        }

        std::vector<Seed> filtered_seeds;
        for (auto it = unitig_to_bucket.begin(); it != unitig_to_bucket.end(); ++it) {
            if (it->second.size() > 1 || seeds[it->second[0]].get_query_view().size() > min_seed_length) {
                filtered_seeds.emplace_back(std::move(seeds[*std::min_element(
                    it.value().begin(), it.value().end(), [&seeds](size_t a, size_t b) {
                        return seeds[a].get_query_view().data() < seeds[b].get_query_view().data();
                    }
                )]));
            }
        }
        common::logger->trace("Done clustering: new seed count: {}", filtered_seeds.size());
        return filtered_seeds;
    }

  private:
    const DeBruijnGraph &graph_;
    sdsl::int_vector<> unitigs_;
    static constexpr auto kUnitigsExtension = ".unitigs";
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __MER_DISTANCES_HPP__
