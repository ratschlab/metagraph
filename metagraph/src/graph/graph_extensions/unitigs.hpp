#ifndef __UNITIGS_HPP__
#define __UNITIGS_HPP__

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "graph/alignment/dbg_aligner.hpp"

namespace mtg {

namespace cli {

class Config;

}

namespace graph {
namespace align {

class Unitigs : public SequenceGraph::GraphExtension {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef bit_vector_small Indicator;
    typedef annot::CountsVector IDVector;
    typedef annot::matrix::TupleRowDiff<annot::ColumnCoordAnnotator::binary_matrix_type> UnitigCoords;
    typedef Alignment::Tuple::value_type Coord;

    Unitigs(const DBGSuccinct &graph) : graph_(std::shared_ptr<const DeBruijnGraph>{}, &graph) {}

    Unitigs(const cli::Config &config);

    bool load(const std::string &filename_base);

    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    void adjacent_outgoing_nodes(size_t unitig_id,
                                 const std::function<void(size_t)> &callback) const {
        graph_->adjacent_outgoing_nodes(get_unitig(unitig_id).second, [&](node_index next) {
            auto next_unitig_ids = get_unitig_ids({ next });
            callback(next_unitig_ids.size() ? next_unitig_ids[0] : next);
        });
    }

    void adjacent_incoming_nodes(size_t unitig_id,
                                 const std::function<void(size_t)> &callback) const {
        graph_->adjacent_incoming_nodes(get_unitig(unitig_id).second, [&](node_index prev) {
            auto prev_unitig_ids = get_unitig_ids({ prev });
            callback(prev_unitig_ids.size() ? prev_unitig_ids[0] : prev);
        });
    }

    std::pair<node_index, node_index> get_unitig(size_t unitig_id) const;

    std::pair<std::pair<node_index, node_index>, std::pair<Coord, Coord>>
    get_unitig_bounds(size_t unitig_id) const;

    std::vector<size_t> get_unitig_ids(const std::vector<node_index> &nodes) const;

    std::vector<std::pair<size_t, Coord>>
    get_unitig_ids_and_coordinates(const std::vector<node_index> &nodes) const;

    size_t cluster_and_filter_seeds(const IDBGAligner &aligner,
                                    IDBGAligner::BatchSeeders &batch_seeders) const;

    void load_graph(const std::string &fname) {
        graph_ = load_graph_impl(fname);
        unitigs_.set_graph(graph_.get());
    }

  private:
    std::shared_ptr<const DBGSuccinct> graph_;
    UnitigCoords unitigs_;
    IDVector boundaries_;
    Indicator indicator_;
    std::unique_ptr<bit_vector> valid_edges_;
    static constexpr auto kUnitigsExtension = ".unitigs";

    static std::shared_ptr<const DBGSuccinct> load_graph_impl(const std::string &fname);

    std::pair<sdsl::bit_vector, std::vector<uint64_t>>
    nodes_to_rows(const std::vector<node_index> &nodes) const;

    size_t get_unitig_id_offset() const { return graph_->max_index() + 1; }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __MER_DISTANCES_HPP__
