#ifndef __PATH_INDEX__HPP
#define __PATH_INDEX__HPP

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"


namespace mtg::graph {

class IPathIndex : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;
    using Row = annot::binmat::BinaryMatrix::Row;
    using RowTuples = annot::matrix::MultiIntMatrix::RowTuples;

    virtual ~IPathIndex() {}

    virtual std::vector<RowTuples> get_coords(const std::vector<node_index> &nodes) const;

  protected:
    virtual size_t coord_to_path_id(uint64_t coord) const = 0;
    virtual uint64_t path_id_to_coord(size_t path_id) const = 0;

    virtual std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows) const = 0;

    virtual bool has_coord(node_index) const { return true; }

    virtual std::pair<size_t, size_t> get_superbubble_terminus(size_t path_id) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
};

template <class PathStorage = annot::RowDiffCoordAnnotator::binary_matrix_type,
          class PathBoundaries = bit_vector_smart,
          class SuperbubbleIndicator = bit_vector_smart>
class PathIndex : public IPathIndex {
  public:
    PathIndex() {}

    PathIndex(std::shared_ptr<const DBGSuccinct> graph,
              const std::string &graph_fname = "",
              const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences
                  = [](const auto &) {});

    PathIndex(const DBGSuccinct &graph,
              const std::string &graph_fname = "",
              const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences
                  = [](const auto &) {})
          : PathIndex(decltype(dbg_succ_)(decltype(dbg_succ_){}, &graph), graph_fname, generate_sequences) {}

    void set_graph(std::shared_ptr<const DBGSuccinct> graph);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool = true) const {
        return dynamic_cast<const DBGSuccinct*>(&graph) == dbg_succ_.get();
    }

  private:
    std::shared_ptr<const DBGSuccinct> dbg_succ_;
    PathStorage paths_indices_;
    PathBoundaries path_boundaries_;
    SuperbubbleIndicator is_superbubble_start_;
    sdsl::int_vector<> superbubble_termini_;

    virtual size_t coord_to_path_id(uint64_t coord) const override final {
        return path_boundaries_.rank1(coord);
    }

    virtual uint64_t path_id_to_coord(size_t path_id) const override final {
        return path_boundaries_.select1(path_id);
    }

    virtual std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows) const override final {
        return paths_indices_.get_row_tuples(rows);
    }

    virtual bool has_coord(node_index node) const override final;

    virtual std::pair<size_t, size_t> get_superbubble_terminus(size_t path_id) const override final;

    virtual const DeBruijnGraph& get_graph() const override final { return *dbg_succ_; }

    static constexpr auto kPathIndexExtension = ".paths";
};

} // namespace mtg::graph

#endif // __PATH_INDEX_HPP__
