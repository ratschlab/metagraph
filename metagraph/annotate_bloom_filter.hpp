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


class AnnotationCategoryBloom : public MultiColorAnnotation<uint64_t, std::string> {
  public:
    using Index = typename MultiColorAnnotation<uint64_t, Color>::Index;
    using Color = typename MultiColorAnnotation<uint64_t, Color>::Color;
    using Coloring = typename MultiColorAnnotation<uint64_t, Color>::Coloring;

    template <typename... Args>
    AnnotationCategoryBloom(const DBG_succ &graph, Args&& ...args)
          : graph_(graph), annotator_(graph_, args...) {}

    void set_coloring(Index i, const Coloring &coloring);
    Coloring get_coloring(Index i) const;

    void add_color(Index i, const Color &color);
    void add_colors(Index i, const Coloring &coloring);
    void add_colors(const std::vector<Index> &indices, const Coloring &coloring);
    void add_colors(const std::string &sequence,
                    const Coloring &colors,
                    size_t num_elements = 0);
    void add_color(const std::string &sequence,
                   const std::string &color,
                   size_t num_elements = 0);

    bool has_color(Index i, const Color &color) const;
    bool has_colors(Index i, const Coloring &coloring) const;

    bool load(const std::vector<std::string> &filenames);
    void serialize(const std::string &filename) const;

    // Get colors that occur at least in |discovery_ratio| colorings.
    // If |discovery_ratio| = 0, return the union of colorings.
    Coloring aggregate_colors(const std::vector<Index> &indices,
                              double discovery_ratio = 1) const;

    // Count all colors collected from extracted colorings
    // and return top |num_top| with the counts computed.
    std::vector<std::pair<Color, size_t>>
    get_most_frequent_colors(const std::vector<Index> &indices,
                             size_t num_top = static_cast<size_t>(-1)) const;

  private:
    DBGSuccAnnotWrapper graph_;
    hash_annotate::BloomAnnotator annotator_;

    std::vector<std::string> column_to_label_;
    std::unordered_map<std::string, uint32_t> label_to_column_;
};

} // namespace annotate

#endif // __ANNOTATE_BLOOM_FILTER_HPP__
