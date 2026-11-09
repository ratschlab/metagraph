#ifndef __TRAVERSAL_TYPES_HPP__
#define __TRAVERSAL_TYPES_HPP__

#include <cstdint>
#include <string>

#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using node_index = SequenceGraph::node_index;
using Column = annot::matrix::BinaryMatrix::Column;
using Row = annot::matrix::BinaryMatrix::Row;
// index of a label in a per-seed label dictionary
using LabelId = uint32_t;

static constexpr node_index npos = SequenceGraph::npos;

// How oriented traversal nodes relate to annotation rows.
enum class Regime {
    BASIC,      // annotation key == node
    PRIMARY,    // CanonicalDBG wrapper over a primary graph: key = get_base_node(node)
    CANONICAL,  // native graph storing both strands: key = the canonical k-mer's node
};
const char* to_string(Regime regime);

// What a label refers to.
enum class LabelKind {
    COLUMN,     // an annotation column (e.g. a sample, file or taxid)
    HEADER,     // one sequence (FASTA header) inside a column, via CoordToHeader
};
const char* to_string(LabelKind kind);

struct LabelRef {
    LabelKind kind = LabelKind::COLUMN;
    Column column = 0;
    uint64_t seq_id = 0;  // only meaningful for HEADER
    std::string name;     // column name or FASTA header
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_TYPES_HPP__
