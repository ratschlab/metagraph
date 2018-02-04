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
    virtual const LabelType& get(Index i) const = 0;
    virtual void set(Index i, const LabelType &label) = 0;

    virtual bool exists(Index i, const LabelType &label) const = 0;

    virtual bool load(const std::string &filename) = 0;
    virtual void serialize(const std::string &filename) const = 0;
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


template <typename LabelType>
class ColorCompressed : public AnnotationCategory<LabelType> {
  public:
    ColorCompressed(const std::vector<ColorCompressed<LabelType>> &categories,
                    const std::vector<std::set<size_t>> &merge_plan);

    virtual const LabelType& get(Index i) const;

    virtual bool exists(Index i, const LabelType &label) const;

    void set(Index i, const LabelType &label);

    virtual bool load(const std::string &filename);
    virtual void serialize(const std::string &filename) const;

    void flush();

  private:
    uint64_t graph_size_;

    std::unordered_map<LabelType, uint64_t> label_to_id_;
    std::vector<LabelType> id_to_label_;
    std::vector<sdsl::sd_vector<>*> bitmatrix_;
    sdsl::bit_vector* annotation_curr_;
    LabelType label_curr_;

    sdsl::bit_vector* inflate_column(const uint64_t id) const;

};


/*
class EdgeWiseMatrix : public UncompressedMatrix {
  public:
};
*/


/*
class WaveletTrie : parent GenomeAnnotation {
  public:
};


class EdgeCompressed : parent GenomeAnnotation {
  public:
};


class ColorBloomFilter : parent GenomeAnnotation {
  public:
};



*/


class DBGSuccAnnotWrapper : public DeBruijnGraphWrapper {
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

    edge_index next_edge(edge_index i, char edge_label) const {
        assert(i >= first_edge() && i <= last_edge());
        return graph_.traverse(i, edge_label);
    }

    edge_index prev_edge(edge_index i) const {
        assert(i >= first_edge() && i <= last_edge());
        assert(has_the_only_incoming_edge(i));
        return graph_.get_minus_k_value(i, 0).second;
    }

  private:
    const DBG_succ &graph_;
};


//TODO: remove this

    // sdsl::bit_vector* inflate_annotation(DBG_succ *G, uint64_t id);

    // void annotate_seq(DBG_succ *G, Config *config, kstring_t &seq, kstring_t &label,
    //                   uint64_t start = 0, uint64_t end = 0,
    //                   pthread_mutex_t *anno_mutex = NULL);

    // // get_annotation(const DBG_succ *G, const std::vector<uint64_t> &node_indices);
    // std::vector<uint32_t> classify_path(DBG_succ *G, std::vector<uint64_t> node_indices);

    // // get_annotation
    // std::set<uint32_t> classify_read(DBG_succ *G, kstring_t &read, uint64_t max_distance);

    // // print
    // void annotationToScreen(DBG_succ *G);

// // write annotation to disk
// void DBG_succ::annotationToFile(const std::string &filename) {
//     std::ofstream outstream(filename);
//     libmaus2::util::NumberSerialisation::serialiseNumber(outstream, annotation_full.size());
//     for (size_t i = 0; i < annotation_full.size(); ++i) {
//         annotation_full.at(i)->serialize(outstream);
//     }
// }

// // read annotation from disk
// void DBG_succ::annotationFromFile(const std::string &filename) {
//     // generate annotation object
//     // populate it with existing annotation if available
//     std::ifstream instream(filename);
//     if (instream.good()) {
//         //if (config->verbose)
//         //    std::cerr << "get annotation from disk" << std::endl;
//         size_t anno_size = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
//         for (size_t i = 0; i < anno_size; ++i) {
//             //annotation_full.push_back(new sdsl::rrr_vector<63>());
//             annotation_full.push_back(new sdsl::sd_vector<>());
//             annotation_full.back()->load(instream);
//         }
//     }
// }

} // namespace annotate

#endif // __ANNOTATE_HPP__
