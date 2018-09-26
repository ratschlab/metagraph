#ifndef __ANNOTATE_BLOOM_FILTER_HPP__
#define __ANNOTATE_BLOOM_FILTER_HPP__

#include <vector>
#include <set>
#include <string>
#include <dbg_bloom_annotator.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"
#include "annotate_column_compressed.hpp"


namespace annotate {

class DBGSuccAnnotWrapper : public hash_annotate::DeBruijnGraphWrapper {
  public:
    DBGSuccAnnotWrapper(const DBG_succ &graph) : graph_(graph) {}

    size_t get_k() const { return graph_.get_k(); }

    edge_index first_edge() const { return 1; }
    edge_index last_edge() const { return graph_.num_nodes(); }

    edge_index map_kmer(const std::string &kmer) const {
        assert(kmer.size() == graph_.get_k());
        return graph_.map_to_nodes(kmer)[0];
    }

    size_t get_num_edges() const {
        return graph_.num_nodes();
    }

    // Transform sequence to the same kind as the de bruijn graph stores
    std::string encode_sequence(const std::string &sequence) const {
        auto result = sequence;
        for (char &c : result) {
            c = graph_.decode(graph_.encode(c));
        }
        return result;
    }

    std::string get_node_kmer(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.get_node_str(i);
    }

    char get_edge_label(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.decode(graph_.get_W(i));
    }

    // Check if the source k-mer for this edge has the only outgoing edge
    bool has_the_only_outgoing_edge(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.is_single_outgoing(i);
    }

    bool has_the_only_incoming_edge(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.indegree(i) == 1;
    }

    bool is_dummy_edge(const std::string &kmer) const {
        assert(kmer.length() == graph_.get_k() + 1);
        return kmer.front() == '$' || kmer.back() == '$';
    }

    bool is_dummy_label(char edge_label) const {
        return edge_label == '$';
    }

    edge_index next_edge(edge_index i, char edge_label) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.traverse(i, edge_label);
    }

    edge_index prev_edge(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        //assert(has_the_only_incoming_edge(i));
        return graph_.bwd(i);
    }

  private:
    const DBG_succ &graph_;
};


class AnnotationCategoryBloom : public MultiLabelAnnotation<uint64_t, std::string> {
  public:
    using Index = typename MultiLabelAnnotation<uint64_t, Label>::Index;
    using Label = typename MultiLabelAnnotation<uint64_t, Label>::Label;
    using VLabels = typename MultiLabelAnnotation<uint64_t, Label>::VLabels;

    template <typename... Args>
    AnnotationCategoryBloom(const DBG_succ &graph, Args&& ...args)
          : graph_(graph), annotator_(graph_, args...) {}

    void set_labels(Index i, const VLabels &labels);
    VLabels get_labels(Index i) const;

    void add_label(Index i, const Label &label);
    void add_labels(Index i, const VLabels &labels);
    void add_labels(const std::vector<Index> &indices, const VLabels &labels);
    void add_labels(const std::string &sequence,
                    const VLabels &labels,
                    size_t num_elements = 0);
    void add_label(const std::string &sequence,
                   const std::string &label,
                   size_t num_elements = 0);

    bool has_label(Index i, const Label &label) const;
    bool has_labels(Index i, const VLabels &labels) const;

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    void insert_rows(const std::vector<Index>&) {}

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    VLabels get_labels(const std::vector<Index> &indices,
                       double presence_ratio) const;

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1)) const;

    uint64_t num_objects() const;
    size_t num_labels() const;
    double sparsity() const;

  private:
    DBGSuccAnnotWrapper graph_;
    hash_annotate::BloomAnnotator annotator_;

    std::vector<std::string> column_to_label_;
    std::unordered_map<std::string, uint32_t> label_to_column_;
};

} // namespace annotate

#endif // __ANNOTATE_BLOOM_FILTER_HPP__
