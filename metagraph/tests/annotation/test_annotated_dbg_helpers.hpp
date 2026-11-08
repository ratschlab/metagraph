#ifndef __TEST_ANNOTATED_DBG_HELPERS_HPP__
#define __TEST_ANNOTATED_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <vector>
#include <string>

#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"


namespace mtg {

namespace test {

class MutableAnnotatedDBG : public graph::AnnotatedDBG {
  public:
    using graph::AnnotatedDBG::AnnotatedDBG;
    using graph::AnnotatedSequenceGraph::graph_;
    using graph::AnnotatedSequenceGraph::annotator_;
};

// |anno| with a sequence header index attached, as an annotation loaded with its .seqs file
inline std::unique_ptr<graph::AnnotatedDBG>
with_coord_to_header(std::unique_ptr<graph::AnnotatedDBG> anno,
                     std::vector<std::vector<std::string>> headers,
                     std::vector<std::vector<uint64_t>> num_kmers) {
    auto &mutable_anno = static_cast<MutableAnnotatedDBG&>(*anno);
    auto graph = std::dynamic_pointer_cast<graph::DeBruijnGraph>(mutable_anno.graph_);
    return std::make_unique<graph::AnnotatedDBG>(
            graph, std::move(mutable_anno.annotator_), false,
            std::make_unique<annot::CoordToHeader>(std::move(headers), std::move(num_kmers)));
}

template <class Graph, class Annotation>
std::unique_ptr<graph::AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences = {},
                 const std::vector<std::string> &labels = {},
                 graph::DeBruijnGraph::Mode mode = graph::DeBruijnGraph::BASIC,
                 bool coordinates = false,
                 const std::vector<uint64_t> &coord_starts = {});

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
