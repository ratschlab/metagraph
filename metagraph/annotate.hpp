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
};

*/


//TODO: remove this
namespace annotate {

    sdsl::bit_vector* inflate_annotation(DBG_succ *G, uint64_t id);

    void annotate_seq(DBG_succ *G, kstring_t &seq, kstring_t &label,
                      uint64_t start = 0, uint64_t end = 0,
                      pthread_mutex_t *anno_mutex = NULL);

    // get_annotation(const DBG_succ *G, const std::vector<uint64_t> &node_indices);
    std::vector<uint32_t> classify_path(DBG_succ *G, std::vector<uint64_t> node_indices);

    // get_annotation
    std::set<uint32_t> classify_read(DBG_succ *G, kstring_t &read, uint64_t max_distance);

    // print
    void annotationToScreen(DBG_succ *G);

} // namespace annotate

#endif // __ANNOTATE_HPP__
