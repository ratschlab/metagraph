#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <map>
#include <memory>

#include "annotated_dbg.hpp"
#include "sequence_graph.hpp"
#include "annotate.hpp"
#include "path.hpp"


class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, AnnotatedDBG::Annotator::VLabels> AlignedPath;

    struct DPAlignmentKey {
        // The node in the underlying graph.
        node_index node;
        std::string::const_iterator query_it;

        bool operator< (const DPAlignmentKey &other) const {
            return node == other.node ? (query_it < other.query_it) : (node < other.node);
        }
    };
    struct DPAlignmentValue {
        // Index of the node which corresponds to the path with lowest cost.
        node_index parent;
        float loss;
    };

    DBGAligner(DeBruijnGraph *graph,
               size_t num_top_paths = 10,
               float sw_threshold = 0.1,
               bool verbose = false,
               float insertion_penalty = 3,
               float deletion_penalty = 3);

    DBGAligner() = delete;
    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    AlignedPath align(const std::string& sequence) const;

    // Return the corresponding sequence of a path according to nodes in the graph.
    std::string get_path_sequence(const std::vector<node_index>& path) const;

  private:
    std::shared_ptr<DeBruijnGraph> graph_;
    // Substitution loss for each pair of nucleotides.
    std::map<char, std::map<char, uint16_t>> sub_loss_;
    // Maximum number of paths to explore at the same time.
    size_t num_top_paths_;
    uint64_t sw_threshold_;
    bool verbose_;
    uint16_t insertion_penalty_;
    uint16_t deletion_penalty_;

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

    // Return the loss of substitution. If not in sub_loss_ return a fixed maximized loss value.
    float single_char_loss(char char_in_query, char char_in_graph) const;

    // Compute the edit distance between the query sequence and the aligned path
    // according to loss parameters in this class.
    float whole_path_loss(const AlignedPath& path, std::string::const_iterator begin) const;

    // Compute the distance between the query sequence and the aligned path sequence
    // according to the CSSW library.
    float ssw_loss(const AlignedPath& path, std::string::const_iterator begin) const;
};

#endif // __DBG_ALIGNER_HPP__
