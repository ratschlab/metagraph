#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <map>
#include <memory>

#include "annotated_dbg.hpp"
#include "sequence_graph.hpp"
#include "path.hpp"


class DBGAligner : public AnnotatedDBG {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, Annotator::VLabels> AlignedPath;
    typedef PathCompare<node_index, Annotator::VLabels> AlignedPathCompare;

    DBGAligner(DeBruijnGraph *dbg,
                 Annotator *annotation,
                 size_t num_threads = 0);

    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    AlignedPath align(const std::string& sequence) const;

    // Compute the edit distance between a node in the graph and a kmer in the string
    // according to loss parameters in this class.
    float single_node_loss(const node_index& node, std::string::const_iterator begin) const;

    // Return the corresponding sequence of a path according to nodes in the graph.
    std::string get_path_sequence(const std::vector<node_index>& path) const;

  private:
    // Substitution loss for each pair of nucleotides. Transition and transversion mutations
    // have different loss values.
    std::map<char, std::map<char, int>> sub_loss_;

    // Align part of a sequence to the graph in the case of no exact map
    // based on internal strategy. Calls callback for every possible alternative path.
    void inexact_map(const AlignedPath &path, std::string::const_iterator end,
                     const std::function<void(node_index,
                        std::string::const_iterator)> &callback) const;

    // Strategy: Randomly choose between possible outgoing neighbors.
    void randomly_pick_strategy(std::vector<node_index> out_neighbors,
                                const std::function<void(node_index)> &callback) const;

    // Strategy: Call callback for all edges. This is equivalent to exhaustive search.
    void pick_all_strategy(std::vector<node_index> out_neighbors,
                           const std::function<void(node_index)> &callback) const;
};

#endif // __DBG_ALIGNER_HPP__
