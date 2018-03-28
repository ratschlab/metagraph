#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <vector>
#include <set>
#include <string>

#include "dbg_succinct.hpp"
#include "dbg_bloom_annotator.hpp"


namespace annotate {

typedef uint64_t Index;


// class GenomeAnnotation {
//   public:
//     typedef uint64_t Index;
//     enum Label {
//         OUTGOING_EDGE_LABELS = 0,
//         SEQUENCE_NAME,
//     };

//     virtual const Annotation& get(Index i, Label j) const = 0;
//     virtual void set(Index i, Label j, const Annotation &label) = 0;
//     virtual void merge(const GenomeAnnotation &annotation) = 0;

//     virtual bool load(const std::string &filename) = 0;
//     virtual void serialize(const std::string &filename) const = 0;
// };


template <typename LabelType>
class AnnotationCategory {
  public:
    virtual ~AnnotationCategory() {}
    virtual LabelType get(Index i) const = 0;
    virtual void set_label(Index i, const LabelType &label) = 0;

    virtual bool has_label(Index i, const LabelType &label) const = 0;

    virtual bool load(const std::string &filename) = 0;
    virtual void serialize(const std::string &filename) const = 0;
};


class ColorCompressed : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    ColorCompressed(uint64_t graph_size)
          : graph_size_(graph_size), annotation_curr_(NULL) {}

    // Merge constructor
    ColorCompressed(const std::vector<ColorCompressed> &categories,
                    const std::vector<std::vector<size_t>> &merge_plan);

    ~ColorCompressed() { release(); }

    SetStr get(Index i) const;

    std::vector<uint64_t> get_row(Index i) const;

    void serialize_uncompressed_rows(const std::string &filename) const;

    bool has_label(Index i, const SetStr &label) const;
    bool has_label(Index i, const std::string &label) const;

    void set_label(Index i, const SetStr &label);
    void add_label(Index i, const std::string &label);

    bool load(const std::string &filename);
    void serialize(const std::string &filename) const;

    std::vector<std::string> get_label_names() const { return id_to_label_; }

  private:
    void release();
    void flush();
    sdsl::bit_vector* inflate_column(const uint32_t id) const;

    uint64_t graph_size_;
    std::unordered_map<std::string, uint32_t> label_to_id_;
    std::vector<std::string> id_to_label_;
    std::vector<sdsl::sd_vector<>*> bitmatrix_;

    uint32_t curr_id_;
    sdsl::bit_vector *annotation_curr_;
};


class PreciseColorCompressedAnnotator : public hash_annotate::PreciseAnnotator {
  public:
    PreciseColorCompressedAnnotator(const ColorCompressed &color_compressed)
          : color_compressed_(color_compressed) {}

    std::vector<uint64_t> annotate_edge(hash_annotate::DeBruijnGraphWrapper::edge_index i) const {
        return color_compressed_.get_row(i);
    }

  private:
    const ColorCompressed &color_compressed_;
};


// class ColorWiseMatrix : public UncompressedMatrix {
//   public:
//     // annotation containers
//     std::deque<uint32_t> annotation; // list that associates each node in the graph with an annotation hash
//     std::vector<std::string> id_to_label; // maps the label ID back to the original string
//     std::unordered_map<std::string, uint32_t> label_to_id_map; // maps each label string to an integer ID
//     std::map<uint32_t, uint32_t> annotation_map; // maps the hash of a combination to the position in the combination vector

//     //std::vector<sdsl::rrr_vector<63>* > annotation_full;
//     std::vector<sdsl::sd_vector<>* > annotation_full;
// };


/*
class EdgeWiseMatrix : public UncompressedMatrix {
  public:
};
*/


/*
class WaveletTrie : public GenomeAnnotation {
  public:
};


class EdgeCompressed : public GenomeAnnotation {
  public:
};

*/


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


class AnnotationCategoryBloom;


class AnnotationCategoryHash : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    AnnotationCategoryHash(const DBG_succ &graph);

    SetStr get(Index i) const;

    void set_label(Index i, const SetStr &label);

    void add_label(const std::string &sequence,
                   const std::string &label);

    bool has_label(Index i, const SetStr &label) const;

    // bool load(const std::string &filename);
    // void serialize(const std::string &filename) const;

    void compare_annotations(const AnnotationCategoryBloom &bloom,
                             size_t step = 1) const;

    void compare_annotations(const hash_annotate::BloomAnnotator &bloom_annotator,
                             size_t step) const;

  private:
    DBGSuccAnnotWrapper graph_;
    hash_annotate::PreciseHashAnnotator annotator_;

    std::vector<std::string> column_to_label_;
    std::unordered_map<std::string, size_t> label_to_column_;
};


class AnnotationCategoryBloom : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    template <typename... Args>
    AnnotationCategoryBloom(const DBG_succ &graph, Args& ...args)
          : graph_(graph), annotator_(graph_, args...) {}

    ~AnnotationCategoryBloom() {}

    SetStr get(Index i) const;

    void set_label(Index i, const SetStr &label);

    void add_label(const std::string &sequence,
                   const std::string &label);

    bool has_label(Index i, const SetStr &label) const;

    bool load(const std::string &filename);
    void serialize(const std::string &filename) const;

    void compare_annotations(const AnnotationCategoryHash &exact,
                             size_t step = 1) const {
        exact.compare_annotations(annotator_, step);
    }

  private:
    DBGSuccAnnotWrapper graph_;
    hash_annotate::BloomAnnotator annotator_;

    std::vector<std::string> column_to_label_;
    std::unordered_map<std::string, uint32_t> label_to_column_;
};

} // namespace annotate

#endif // __ANNOTATE_HPP__
