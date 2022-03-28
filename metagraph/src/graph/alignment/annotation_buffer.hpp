#ifndef __ANNOTATION_BUFFER_HPP__
#define __ANNOTATION_BUFFER_HPP__

#include "alignment.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vector_set.hpp"
#include "common/vector_map.hpp"
#include "common/hashers/hash.hpp"

namespace mtg {
namespace graph {

class CanonicalDBG;

namespace align {

// caches queried annotations to speed up next queries (labels with or w/o coordinates)
class AnnotationBuffer {
  public:
    typedef AnnotatedDBG::Annotator Annotator;
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::Tuple Tuple;
    typedef Alignment::Columns Columns;
    typedef Alignment::CoordinateSet CoordinateSet;

    AnnotationBuffer(const DeBruijnGraph &graph,
                     const Annotator &annotator,
                     size_t row_batch_size = std::numeric_limits<size_t>::max());

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

    // get the labels of a node if they have been fetched
    inline const Columns* get_labels(node_index node) const {
        return get_labels_and_coords(node).first;
    }

    const Annotator& get_annotator() const { return annotator_; }

    size_t num_nodes_buffered() const { return node_to_cols_.size(); }
    size_t num_column_sets() const { return column_sets_.size(); }

    // This method lets the caller push additional column sets to dictionary
    // `column_sets_`. These column sets can later be fetched with `get_cached_column_set()`
    template <typename... Args>
    inline size_t cache_column_set(Args&&... args) {
        auto it = column_sets_.emplace(std::forward<Args>(args)...).first;
        assert(std::is_sorted(it->begin(), it->end()));
        return it - column_sets_.begin();
    }

    // Fetch a label set given its index returned by `cache_column_set()`
    inline const Columns& get_cached_column_set(size_t i) const {
        assert(i < column_sets_.size());
        return column_sets_.data()[i];
    }

  private:
    const DeBruijnGraph &graph_;
    const Annotator &annotator_;
    const annot::matrix::MultiIntMatrix *multi_int_;
    const CanonicalDBG *canonical_;

    // keep a unique set of annotation rows
    // the first element is the empty label set
    VectorSet<Columns, utils::VectorHash> column_sets_;
     // map node to index in |column_sets_|
    VectorMap<node_index, size_t> node_to_cols_;
    // coordinate sets for all nodes in |node_to_cols_| in the same order
    std::vector<CoordinateSet> label_coords_;
    // buffer of paths to later querying with fetch_queued_annotations()
    std::vector<std::vector<node_index>> queued_paths_;

    size_t row_batch_size_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ANNOTATION_BUFFER_HPP__
