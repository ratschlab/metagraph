#ifndef __TRAVERSAL_TYPES_HPP__
#define __TRAVERSAL_TYPES_HPP__

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using node_index = SequenceGraph::node_index;
using Column = annot::matrix::BinaryMatrix::Column;
using Row = annot::matrix::BinaryMatrix::Row;
using Coord = uint64_t;
// index of a label in a per-seed label dictionary
using LabelId = uint32_t;

static constexpr node_index npos = SequenceGraph::npos;

// How oriented traversal nodes relate to annotation rows.
enum class Regime {
    BASIC,      // annotation key == node
    PRIMARY,    // CanonicalDBG wrapper over a primary graph: key = get_base_node(node)
    CANONICAL,  // native graph storing both strands: key = map_to_nodes(k-mer)
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

    bool same_target(const LabelRef &other) const {
        return kind == other.kind && column == other.column
            && (kind == LabelKind::COLUMN || seq_id == other.seq_id);
    }
};

// Evidence required for a label to support a k-mer on a path.
enum class Support {
    KMER,   // the label carries every k-mer of the path (label-consistent)
    TRACE,  // consecutive k-mer coordinates in the same indexed sequence (trace-consistent)
};
const char* to_string(Support support);

enum class Arm : uint8_t { LEFT = 0, RIGHT = 1 };
const char* to_string(Arm arm);

// Why a label left a path, or why a path ended.
enum class EndReason : uint8_t {
    // semantic, per label
    DEAD_END = 0,     // no structural successor in this direction
    LABEL_LOST,       // successors exist but none carries the label
    LOSS_BUDGET,      // a label switch was possible but only above the loss budget
    BRANCH,           // the label is on several successors and its branch limit is exhausted
    EDGE_REUSE,       // the label's continuation uses an edge already used on this path
    EDGE_REUSE_RC,    // ... used in the opposite orientation (hairpin / inverted repeat)
    REACHED_SEED,     // ... the edge belongs to the seed (e.g. a closed circle)
    RECORD_END,       // trace mode: the indexed source sequence ends
    MAX_EXTENSION,    // the requested radius was reached (requested domain complete)
    // resource caps: exploration of the requested domain is incomplete
    MAX_STEPS,
    MAX_LIVE_PATHS,
    MAX_PATHS,
    MAX_OUTPUT,
    TIME_BUDGET,
    BEAM_PRUNED,
};
static constexpr size_t kNumEndReasons = static_cast<size_t>(EndReason::BEAM_PRUNED) + 1;
const char* to_string(EndReason reason);
// true if the reason means the requested domain was not explored completely
bool is_resource_stop(EndReason reason);

static constexpr double kInfiniteLoss = std::numeric_limits<double>::infinity();

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_TYPES_HPP__
