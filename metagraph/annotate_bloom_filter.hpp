#ifndef __ANNOTATE_BLOOM_FILTER_HPP__
#define __ANNOTATE_BLOOM_FILTER_HPP__

#include <vector>
#include <set>
#include <string>
#include <annograph/bloom_annotator/dbg_bloom_annotator.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"
#include "annotate_color_compressed.hpp"


namespace annotate {

class DBGSuccAnnotWrapper : public hash_annotate::DeBruijnGraphWrapper {
  public:
    DBGSuccAnnotWrapper(const DBG_succ &graph) : graph_(graph) {}

    size_t get_k() const { return graph_.get_k(); }

    edge_index first_edge() const { return 1; }
    edge_index last_edge() const { return graph_.get_W().size() - 1; }

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
        return graph_.get_last(i) && ((i == 1) || graph_.get_last(i - 1));
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
        return graph_.get_minus_k_value(i, 0).second;
    }

  private:
    const DBG_succ &graph_;
};


class AnnotationCategoryHash;


class AnnotationCategoryBloom : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    template <typename... Args>
    AnnotationCategoryBloom(const DBG_succ &graph, Args&& ...args)
          : graph_(graph), annotator_(graph_, args...) {}

    ~AnnotationCategoryBloom() {}

    SetStr get(Index i) const;

    void set_label(Index i, const SetStr &label);

    void add_label(const std::string &sequence,
                   const std::string &label,
                   size_t num_elements = 0);

    void add_labels(const std::string &sequence,
                    const SetStr &labels,
                    size_t num_elements);

    void add_labels(const std::string &sequence,
                    const SetStr &labels) {
        add_labels(sequence, labels, 0);
    }

    bool has_label(Index i, const SetStr &label) const;

    bool load(const std::string &filename);
    void serialize(const std::string &filename) const;

    std::vector<std::string> get_label_names() const {
        return column_to_label_;
    }

  private:
    DBGSuccAnnotWrapper graph_;
    hash_annotate::BloomAnnotator annotator_;

    std::vector<std::string> column_to_label_;
    std::unordered_map<std::string, uint32_t> label_to_column_;
};

} // namespace annotate

#endif // __ANNOTATE_BLOOM_FILTER_HPP__
