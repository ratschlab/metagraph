#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <cassert>
#include <memory>
#include <mutex>

#include <sdsl/int_vector.hpp>

#include "representation/base/sequence_graph.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace graph {

class AnnotatedSequenceGraph {
  public:
    typedef std::string Label;
    typedef annot::MultiLabelEncoded<Label> Annotator;
    using node_index = SequenceGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graph,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast = false);

    virtual ~AnnotatedSequenceGraph() {}

    virtual std::vector<Label> get_labels(node_index index) const;

    virtual bool has_label(node_index index, const Label &label) const;

    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequence(std::string_view sequence,
                                   const std::vector<Label> &labels);
    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequences(
        const std::vector<std::pair<std::string, std::vector<Label>>> &data);

    virtual void call_annotated_nodes(const Label &label,
                                      std::function<void(node_index)> callback) const;

    virtual bool label_exists(const Label &label) const;

    virtual bool check_compatibility() const;

    virtual const SequenceGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const SequenceGraph> get_graph_ptr() const { return graph_; }

    virtual const Annotator& get_annotator() const { return *annotator_; }

    static row_index graph_to_anno_index(node_index kmer_index) {
        assert(kmer_index);
        return kmer_index - 1;
    }
    static node_index anno_to_graph_index(row_index anno_index) {
        return anno_index + 1;
    }

  protected:
    std::shared_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    std::mutex mutex_;
    bool force_fast_;
};


class AnnotatedDBG : public AnnotatedSequenceGraph {
  public:
    AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 bool force_fast = false);

    using AnnotatedSequenceGraph::get_labels;

    const DeBruijnGraph& get_graph() const { return dbg_; }

    // add k-mer counts to the annotation, thread-safe for concurrent calls
    void add_kmer_counts(std::string_view sequence,
                         const std::vector<Label> &labels,
                         std::vector<uint64_t>&& kmer_counts);

    // add k-mer coordinates to the annotation, the binary annotation must exist
    void add_kmer_coord(std::string_view sequence,
                        const std::vector<Label> &labels,
                        uint64_t start);

    // add k-mer coordinates to the annotation, the binary annotation must exist
    void add_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data);

    // annotate k-mer and their coordinates (combines annotate_sequences and add_kmer_coords)
    void annotate_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data);

    /*********************** Special queries **********************/

    // return labels that occur at least in |discovery_fraction| k-mers
    // but skip the sequence if fewer than |presence_fraction| k-mers are matched against the graph
    std::vector<Label> get_labels(std::string_view sequence,
                                  double discovery_fraction = 0.0,
                                  double presence_fraction = 0.0) const;

    std::vector<Label> get_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                                  size_t min_count) const;

    // Return top |num_top_labels| labels with their counts.
    // The returned counts are weighted by the annotated relation counts if
    // |with_kmer_counts| is true.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(std::string_view sequence,
                   size_t num_top_labels,
                   double discovery_fraction = 0.0,
                   double presence_fraction = 0.0,
                   bool with_kmer_counts = false) const;

    // The returned counts are weighted by the annotated relation counts if
    // |with_kmer_counts| is true.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                   size_t num_top_labels,
                   size_t min_count = 0,
                   bool with_kmer_counts = false) const;

    // returns tuples (label, num_kmer_matches, kmer_abundance_quantiles)
    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>
    get_label_count_quantiles(std::string_view sequence,
                              size_t num_top_labels,
                              double discovery_fraction,
                              double presence_fraction,
                              const std::vector<double> &count_quantiles) const;

    // returns tuples (label, num_kmer_matches, kmer_abundances)
    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>
    get_kmer_counts(std::string_view sequence,
                    size_t num_top_labels,
                    double discovery_fraction,
                    double presence_fraction) const;

    // returns tuples (label, num_kmer_matches, kmer_abundances)
    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>
    get_kmer_counts(const std::vector<node_index> &nodes,
                    size_t num_top_labels,
                    double discovery_fraction,
                    double presence_fraction) const;

    // returns tuples (label, num_kmer_matches, kmer_coordinates)
    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>
    get_kmer_coordinates(std::string_view sequence,
                         size_t num_top_labels,
                         double discovery_fraction,
                         double presence_fraction) const;

    // returns tuples (label, num_kmer_matches, kmer_coordinates)
    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>
    get_kmer_coordinates(const std::vector<node_index> &nodes,
                         size_t num_top_labels,
                         double discovery_fraction,
                         double presence_fraction) const;

    std::vector<std::pair<Label, sdsl::bit_vector>>
    get_top_label_signatures(std::string_view sequence,
                             size_t num_top_labels,
                             double discovery_fraction = 0.0,
                             double presence_fraction = 0.0) const;

    int32_t score_kmer_presence_mask(const sdsl::bit_vector &kmer_presence_mask,
                                     int32_t match_score = 1,
                                     int32_t mismatch_score = 2) const;

  private:
    DeBruijnGraph &dbg_;
};

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_DBG_HPP__
