#ifndef __PATH_INDEX__HPP
#define __PATH_INDEX__HPP

#include <sdsl/dac_vector.hpp>
#include <cache.hpp>
#include <lru_cache_policy.hpp>
#include <tsl/hopscotch_set.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg::graph {

class ColumnPathIndex : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;
    using Row = annot::binmat::BinaryMatrix::Row;
    using RowTuples = annot::matrix::MultiIntMatrix::RowTuples;
    using Label = std::string;
    using SequenceCallback = std::function<void(const std::string&)>;
    using SequenceGenerator = std::function<void(const SequenceCallback&)>;
    using FastaCallback = std::function<void(const std::string &, const SequenceGenerator&)>;
    using FastaGenerator = std::function<void(const FastaCallback&)>;

    ColumnPathIndex(std::shared_ptr<const DBGSuccinct> graph,
                    const annot::LabelEncoder<Label> &label_encoder,
                    const FastaGenerator &fasta_generator,
                    const std::string &graph_fname = "");

    ColumnPathIndex(const DBGSuccinct &graph,
                    const annot::LabelEncoder<Label> &label_encoder,
                    const FastaGenerator &fasta_generator,
                    const std::string &graph_fname = "")
          : ColumnPathIndex(decltype(dbg_succ_)(decltype(dbg_succ_){}, &graph),
                            label_encoder,
                            fasta_generator,
                            graph_fname) {}

    std::vector<RowTuples> get_coords(const std::vector<node_index> &nodes) const;

    void adjacent_outgoing_unitigs(size_t path_id,
                                   const std::function<void(size_t)> &callback) const;

    void set_graph(std::shared_ptr<const DBGSuccinct> graph);

    bool load(const std::string &) { return true; }
    void serialize(const std::string &) const {}

    bool is_compatible(const SequenceGraph &graph, bool = true) const {
        return dynamic_cast<const DBGSuccinct*>(&graph) == dbg_succ_.get();
    }

  private:
    std::shared_ptr<const DBGSuccinct> dbg_succ_;
};

class IPathIndex : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;
    using Row = annot::binmat::BinaryMatrix::Row;
    using RowTuples = annot::matrix::MultiIntMatrix::RowTuples;
    using SuperbubbleStorage = sdsl::dac_vector_dp<>;
    using SuperbubbleInfo = std::tuple<size_t,
                                       SuperbubbleStorage::const_iterator,
                                       SuperbubbleStorage::const_iterator>;

    virtual ~IPathIndex() {}

    virtual std::vector<RowTuples> get_coords(const std::vector<node_index> &nodes) const;

    virtual SuperbubbleInfo
    get_superbubble_terminus(size_t path_id) const = 0;

    virtual SuperbubbleInfo
    get_superbubble_and_dist(size_t path_id) const = 0;

    virtual SuperbubbleInfo get_superbubble_chain(size_t path_id) const = 0;

    virtual size_t coord_to_path_id(uint64_t coord) const = 0;
    virtual size_t coord_to_read_id(uint64_t coord) const = 0;
    virtual uint64_t path_id_to_coord(size_t path_id) const = 0;
    virtual uint64_t read_id_to_coord(size_t read_id) const = 0;
    virtual bool can_reach_superbubble_terminus(size_t path_id) const = 0;

    virtual void adjacent_outgoing_unitigs(size_t path_id,
                                           const std::function<void(size_t)> &callback) const = 0;

    virtual bool is_unitig(size_t path_id) const = 0;

    size_t path_length(size_t path_id) const {
        assert(is_unitig(path_id));
        return path_id_to_coord(path_id + 1) - path_id_to_coord(path_id);
    }

    size_t read_length(size_t read_id) const {
        assert(!is_unitig(read_id));
        return read_id_to_coord(read_id + 1) - read_id_to_coord(read_id);
    }

    void call_dists(size_t path_id_1,
                    size_t path_id_2,
                    const std::function<void(size_t)> &callback,
                    size_t max_dist = std::numeric_limits<size_t>::max(),
                    size_t max_search_depth = 2) const;

  protected:
    virtual std::vector<RowTuples> get_path_row_tuples(const std::vector<Row> &rows) const = 0;
    virtual std::vector<RowTuples> get_read_row_tuples(const std::vector<Row> &rows) const = 0;

    virtual bool has_coord(node_index) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
};

template <class PathStorage = annot::RowDiffCoordAnnotator::binary_matrix_type,
          class PathBoundaries = bit_vector_smart,
          class SuperbubbleIndicator = bit_vector_smart,
          class SuperbubbleStorage = IPathIndex::SuperbubbleStorage>
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

    virtual SuperbubbleInfo
    get_superbubble_terminus(size_t path_id) const override final;

    virtual SuperbubbleInfo
    get_superbubble_and_dist(size_t path_id) const override final;

    virtual bool can_reach_superbubble_terminus(size_t path_id) const override final;

    virtual SuperbubbleInfo get_superbubble_chain(size_t path_id) const override final;

    virtual size_t coord_to_path_id(uint64_t coord) const override final {
        return path_boundaries_.rank1(coord);
    }

    virtual size_t coord_to_read_id(uint64_t coord) const override final {
        return num_unitigs_ + read_boundaries_.rank1(coord);
    }

    virtual uint64_t path_id_to_coord(size_t path_id) const override final {
        return path_boundaries_.select1(path_id);
    }

    virtual uint64_t read_id_to_coord(size_t read_id) const override final {
        return read_boundaries_.select1(read_id);
    }

    virtual void adjacent_outgoing_unitigs(size_t path_id,
                                           const std::function<void(size_t)> &callback) const override final;

    virtual bool is_unitig(size_t path_id) const override final {
        return path_id && path_id <= num_unitigs_;
    }

  private:
    std::shared_ptr<const DBGSuccinct> dbg_succ_;
    size_t num_unitigs_;
    size_t num_superbubbles_;
    PathBoundaries dummy_indicator_;
    PathStorage paths_indices_;
    PathBoundaries path_boundaries_;

    PathStorage read_indices_;
    PathBoundaries read_boundaries_;

    SuperbubbleStorage unitig_backs_;
    SuperbubbleStorage unitig_fronts_;

    SuperbubbleStorage superbubble_sources_;
    SuperbubbleIndicator superbubble_sources_b_;

    SuperbubbleStorage superbubble_termini_;
    SuperbubbleIndicator superbubble_termini_b_;

    SuperbubbleIndicator can_reach_terminus_;

    SuperbubbleStorage unitig_chain_;
    SuperbubbleStorage unitig_chain_sizes_;
    SuperbubbleIndicator unitig_chain_b_;

    virtual std::vector<RowTuples> get_path_row_tuples(const std::vector<Row> &rows) const override final {
        return paths_indices_.get_row_tuples(rows);
    }
    virtual std::vector<RowTuples> get_read_row_tuples(const std::vector<Row> &rows) const override final {
        return read_indices_.num_rows()
            ? read_indices_.get_row_tuples(rows)
            : std::vector<RowTuples>(rows.size(), RowTuples{});
    }

    virtual bool has_coord(node_index node) const override final {
        assert(node < dummy_indicator_.size());
        return node != DeBruijnGraph::npos && !dummy_indicator_[node];
    }

    virtual const DeBruijnGraph& get_graph() const override final { return *dbg_succ_; }

    static constexpr auto kPathIndexExtension = ".paths";
};

} // namespace mtg::graph

#endif // __PATH_INDEX_HPP__
