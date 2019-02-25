#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <memory>

#include "annotated_dbg.hpp"
#include "sequence_graph.hpp"


class DBGAligner : public AnnotatedDBG {
  public:
    typedef DeBruijnGraph::node_index node_index;

    DBGAligner(DeBruijnGraph *dbg,
                 Annotator *annotation,
                 size_t num_threads = 0);

    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    std::vector<node_index> align(const std::string& sequence) const;

    // Return the corresponding sequence of a path according to nodes in the graph.
    std::string get_path_sequence(const std::vector<node_index>& path) const;

  private:
    // Align part of a sequence to the graph in the case of no exact map
    // based on internal strategy. Terminates at the end of sequence or
    // as soon as there is an exact match to a node in graph.
    void inexact_map(std::string::const_iterator begin, std::string::const_iterator end,
                     const node_index &last_mapped_node, const Annotator::VLabels &label_set,
                     const std::function<void(node_index)> &callback,
                     const std::function<bool()> &terminate = []() { return false; }) const;
};

#endif // __DBG_ALIGNER_HPP__
