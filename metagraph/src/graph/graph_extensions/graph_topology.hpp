#ifndef __GRAPH_TOPOLOGY_HPP__
#define __GRAPH_TOPOLOGY_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/logger.hpp"


namespace mtg::graph {

class GraphTopology : public SequenceGraph::GraphExtension {
  public:
    using Coords = Vector<std::tuple<annot::matrix::BinaryMatrix::Column, size_t, size_t, uint64_t>>;

    GraphTopology(std::unique_ptr<annot::matrix::MultiIntMatrix>&& coords,
                  std::unique_ptr<annot::matrix::BinaryMatrix>&& unitigs,
                  std::unique_ptr<annot::matrix::BinaryMatrix>&& clusters)
          : coords_(std::move(coords)),
            unitig_indicator_(std::move(unitigs)),
            cluster_indicator_(std::move(clusters)) {}

    std::vector<Coords> get_coords(const std::vector<annot::matrix::BinaryMatrix::Row> &rows);

    bool load(const std::string &filename_base) {
        {
            std::ifstream fin(filename_base + kCoordExtension, std::ios::binary);
            if (!coords_->load(fin)) {
                common::logger->error("Failed to load coords from {}", filename_base + kCoordExtension);
                return false;
            }
        }

        {
            std::ifstream fin(filename_base + kUnitigExtension, std::ios::binary);
            if (!unitig_indicator_->load(fin)) {
                common::logger->error("Failed to load unitig indicator from {}", filename_base + kUnitigExtension);
                return false;
            }
        }

        {
            std::ifstream fin(filename_base + kClusterExtension, std::ios::binary);
            if (!cluster_indicator_->load(fin)) {
                common::logger->error("Failed to load cluster indicator from {}", filename_base + kClusterExtension);
                return false;
            }
        }

        return true;
    }

    void serialize(const std::string &filename_base) const {
        {
            std::ofstream fout(filename_base + kCoordExtension, std::ios::binary);
            coords_->serialize(fout);
        }
        {
            std::ofstream fout(filename_base + kUnitigExtension, std::ios::binary);
            unitig_indicator_->serialize(fout);
        }
        {
            std::ofstream fout(filename_base + kClusterExtension, std::ios::binary);
            cluster_indicator_->serialize(fout);
        }
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

  private:
    std::unique_ptr<annot::matrix::MultiIntMatrix> coords_;
    std::unique_ptr<annot::matrix::BinaryMatrix> unitig_indicator_;
    std::unique_ptr<annot::matrix::BinaryMatrix> cluster_indicator_;

    static constexpr auto kCoordExtension = ".coord.annodbg";
    static constexpr auto kUnitigExtension = ".unitig.annodbg";
    static constexpr auto kClusterExtension = ".cluster.annodbg";
};

} // namespace mtg::graph

#endif // __GRAPH_TOPOLOGY_HPP__
