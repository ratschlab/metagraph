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
    using Label = AnnotatedDBG::Label;
    using ChainInfo = std::tuple<size_t, size_t, size_t>;
    using CoordInfo = std::tuple<int64_t, Vector<int64_t>, Vector<int64_t>>;
    using InfoPair = std::pair<ChainInfo, CoordInfo>;
    using NodesInfo = Vector<InfoPair>;
    using LabeledNodesInfo = std::pair<Label, NodesInfo>;

    using Labels = std::vector<Label>;
    using LabeledSeqCallback = std::function<void(std::string, size_t, size_t, size_t, std::string, const Labels&)>;
    using LabeledSeqGenerator = std::function<void(const LabeledSeqCallback&)>;

    ColumnPathIndex(std::shared_ptr<DeBruijnGraph> graph,
                    const std::vector<std::string> &prefixes);

    static void annotate_columns(std::shared_ptr<DeBruijnGraph> graph,
                                 const std::string &out_prefix,
                                 const LabeledSeqGenerator &generator,
                                 size_t num_columns_cached = 10,
                                 const std::string &tmp_dir = "",
                                 double memory_available_gb = 1.0,
                                 size_t max_chunks_open = 2000);

    ColumnPathIndex(std::shared_ptr<const AnnotatedDBG> anno_graph,
                    std::shared_ptr<const AnnotatedDBG::Annotator> topo_annotator);

    ColumnPathIndex(const AnnotatedDBG &anno_graph,
                    const AnnotatedDBG::Annotator &topo_annotator)
          : anno_graph_(std::shared_ptr<const AnnotatedDBG>{}, &anno_graph),
            topo_annotator_(std::shared_ptr<const AnnotatedDBG::Annotator>{}, &topo_annotator) {}

    InfoPair get_chain_info(const Label &label, node_index node) const;
    std::vector<LabeledNodesInfo> get_chain_info(const std::vector<node_index> &nodes) const;

    void call_distances(const Label &label,
                        const InfoPair &info_a,
                        const InfoPair &info_b,
                        const std::function<void(size_t)> &callback,
                        int64_t max_distance = std::numeric_limits<int64_t>::max(),
                        size_t max_steps = std::numeric_limits<int64_t>::max()) const;

    bool load(const std::string &) { throw std::runtime_error("Load not implemented for ColumnPathIndex"); }

    void serialize(const std::string &filename) const {
        anno_graph_->get_annotator().serialize(filename + ".global");
        topo_annotator_->serialize(filename + ".topo");
    }

    bool is_compatible(const SequenceGraph &graph, bool = true) const {
        return anno_graph_->check_compatibility() && &graph == &anno_graph_->get_graph();
    }

    const AnnotatedDBG& get_anno_graph() const { return *anno_graph_; }
    const AnnotatedDBG::Annotator& get_topology_annotator() const { return *topo_annotator_; }

    static const std::string UNITIG_FRONT_TAG;
    static const std::string UNITIG_BACK_TAG;
    static const std::string SUPERBUBBLE_TAG;
    static const std::string CHAIN_TAG;


  private:
    std::shared_ptr<const AnnotatedDBG> anno_graph_;
    std::shared_ptr<const AnnotatedDBG::Annotator> topo_annotator_;

    ColumnPathIndex(std::pair<std::shared_ptr<const AnnotatedDBG>,
                              std::shared_ptr<const AnnotatedDBG::Annotator>> annotators)
          : ColumnPathIndex(annotators.first, annotators.second) {}

    const annot::matrix::IntMatrix& get_topo_matrix() const;
    const annot::matrix::MultiIntMatrix& get_coord_matrix() const;

    int64_t get_global_coord(const Label &label, size_t unitig_id) const;
    node_index get_unitig_back(const Label &label, size_t unitig_id) const;

    void adjacent_outgoing_unitigs(const Label &label, size_t unitig_id, const std::function<void(size_t, int64_t)> &callback) const;
};

} // namespace mtg::graph

#endif // __PATH_INDEX_HPP__
