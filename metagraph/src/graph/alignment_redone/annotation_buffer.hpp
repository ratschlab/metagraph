#pragma once

#include "graph/annotated_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vector_set.hpp"
#include "common/vector_map.hpp"
#include "common/hashers/hash.hpp"

namespace mtg::graph {

class CanonicalDBG;

namespace align_redone {

// caches queried annotations to speed up next queries (labels with or w/o coordinates)
class AnnotationBuffer {
  public:
    using label_class_t = size_t;
    using coord_t = int64_t;

    using Annotator = AnnotatedDBG::Annotator;
    using node_index = DeBruijnGraph::node_index;

    using Column = annot::matrix::BinaryMatrix::Column;
    using Tuple = SmallVector<coord_t>;
    using Columns = annot::matrix::BinaryMatrix::SetBitPositions;
    using CoordinateSet = Vector<Tuple>;

    // dummy index for an unfetched annotations
    static constexpr label_class_t nannot = std::numeric_limits<label_class_t>::max();
    static constexpr Column ncolumn = std::numeric_limits<Column>::max();
    static constexpr coord_t ncoord = std::numeric_limits<coord_t>::max();

    AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator);

    void queue_path(const std::vector<node_index> &path) {
        queued_paths_.push_back(path);
    }

    void queue_path(std::vector<node_index>&& path) {
        queued_paths_.push_back(std::move(path));
    }

    // fetch annotations for the queued nodes from the buffer and reset the buffer
    void fetch_queued_annotations();

    bool has_coordinates() const { return multi_int_; }

    // Get the annotations and coordinates of a node if they have been fetched.
    // The returned pointers are valid until next fetch_queued_annotations().
    std::pair<const Columns*, const CoordinateSet*>
    get_labels_and_coords(node_index node) const;

    std::pair<label_class_t, const CoordinateSet*>
    get_labels_id_and_coords(node_index node) const;

    // get the labels of a node if they have been fetched
    inline const Columns* get_labels(node_index node) const {
        return get_labels_and_coords(node).first;
    }

    // get the labels of a node if they have been fetched
    inline label_class_t get_labels_id(node_index node) const {
        return get_labels_id_and_coords(node).first;
    }

    const Annotator& get_annotator() const { return annotator_; }

    size_t num_nodes_buffered() const { return node_to_cols_.size(); }
    size_t num_column_sets() const { return column_sets_.size(); }

    // This method lets the caller push additional column sets to dictionary
    // `column_sets_`. These column sets can later be fetched with `get_cached_column_set()`
    template <typename... Args>
    inline label_class_t cache_column_set(Args&&... args) {
        auto it = column_sets_.emplace(std::forward<Args>(args)...).first;
        assert(std::is_sorted(it->begin(), it->end()));
        return it - column_sets_.begin();
    }

    inline label_class_t cache_column(Column column) {
        return cache_column_set(Columns(1, column));
    }

    // Fetch a label set given its index returned by `cache_column_set()`
    inline const Columns& get_cached_column_set(label_class_t i) const {
        assert(i < column_sets_.size());
        return column_sets_.data()[i];
    }

    std::string generate_column_set_str(label_class_t i, size_t spelling_size) const;

    CoordinateSet get_label_path_coords(const std::vector<node_index> &path,
                                        const std::vector<label_class_t> &label_path) const;

  private:
    const DeBruijnGraph &graph_;
    const Annotator &annotator_;
    const annot::matrix::MultiIntMatrix *multi_int_;
    const CanonicalDBG *canonical_;

    // keep a unique set of annotation rows
    // the first element is the empty label set
    VectorSet<Columns, utils::VectorHash> column_sets_;
     // map node to index in |column_sets_|
    VectorMap<node_index, label_class_t> node_to_cols_;
    // coordinate sets for all nodes in |node_to_cols_| in the same order
    std::vector<CoordinateSet> label_coords_;
    // buffer of paths to later querying with fetch_queued_annotations()
    std::vector<std::vector<node_index>> queued_paths_;
};

} // namespace align
} // namespace mtg::graph
