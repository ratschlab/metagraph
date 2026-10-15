#ifndef __CONTEXT_CURSOR_HPP__
#define __CONTEXT_CURSOR_HPP__

#include <array>
#include <string>
#include <vector>

#include "graph/annotated_dbg.hpp"

namespace mtg::graph { class NodeFirstCache; }

namespace mtg::graph::traversal {

// An exact substring context over indexed source records. Both directions share
// one occurrence set; choosing a flank also constrains the opposite flank.
// Requires an immutable BASIC DNA graph, coordinates, and record boundaries.
// The graph and annotation must outlive the cursor. A cursor is not thread-safe.
class ContextCursor {
  public:
    enum class Side { LEFT, RIGHT };
    struct Occurrence {
        uint64_t column;
        size_t record;
        uint64_t offset; // zero-based start of the current sequence in the record
    };
    struct Choice {
        char left;  // '\0': record end; '?': continuation absent from the index
        char right;
        size_t count;
    };

    ContextCursor(const AnnotatedDBG &graph, std::string seed,
                  const std::vector<std::string> &columns = {},
                  size_t max_occurrences = 100000);
    // Branch exploration can copy a state without decoding its seed again.
    ContextCursor(const ContextCursor &other);
    ContextCursor(ContextCursor &&other) noexcept;
    ~ContextCursor();

    const std::string &sequence() const { return sequence_; }
    std::vector<Occurrence> occurrences() const;
    std::vector<Choice> choices();

    // Restrict to occurrences supporting this character, then extend by one.
    // Returns false, without changing the cursor, if the choice is absent.
    bool extend(Side side, char base);
    // Extend each flank while every occurrence supports the same character.
    // Record ends stop closure. Returns (left extension, right extension).
    std::pair<size_t, size_t> close(size_t max_bases_per_side);

  private:
    using Node = DeBruijnGraph::node_index;
    struct Hit {
        Occurrence occurrence;
        uint64_t global;
        uint64_t record_kmers;
    };
    void find_boundaries();

    const AnnotatedDBG &annotated_;
    std::string sequence_;
    Node first_ = 0, last_ = 0;
    std::vector<Hit> hits_;
    std::array<bool, 2> boundaries_ready_{};
    std::unique_ptr<NodeFirstCache> first_cache_;
    std::vector<std::array<char, 2>> boundaries_;
    std::array<std::array<Node, 4>, 2> neighbors_{};
};

} // namespace mtg::graph::traversal
#endif
