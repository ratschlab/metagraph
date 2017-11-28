#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <vector>
#include <set>
#include <string>

#include "dbg_succinct.hpp"


/*
class GenomeAnnotation {
  public:
      get_annotation = 0;
      get_annotation
      get_annotated_edges(G, annotation)

      serialize
      load

      print
      std::ostream& operator<<();
};


class WaveletTrie : parent GenomeAnnotation {
  public:
};


class ColorCompressed : parent GenomeAnnotation {
  public:
};


class EdgeCompressed : parent GenomeAnnotation {
  public:
};


class ColorBloomFilter : parent GenomeAnnotation {
  public:
};


class UncompressedMatrix : parent GenomeAnnotation {
  public:
};


class EdgeWiseMatrix : parent UncompressedMatrix {
  public:
};


class ColorWiseMatrix : parent UncompressedMatrix {
  public:
    // annotation containers
    std::deque<uint32_t> annotation; // list that associates each node in the graph with an annotation hash
    std::vector<std::string> id_to_label; // maps the label ID back to the original string
    std::unordered_map<std::string, uint32_t> label_to_id_map; // maps each label string to an integer ID
    std::map<uint32_t, uint32_t> annotation_map; // maps the hash of a combination to the position in the combination vector

    //std::vector<sdsl::rrr_vector<63>* > annotation_full;
    std::vector<sdsl::sd_vector<>* > annotation_full;
};

*/


//TODO: remove this
namespace annotate {

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
