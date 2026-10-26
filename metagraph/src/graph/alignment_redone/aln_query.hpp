#pragma once

#include <string>
#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg::graph::align_redone {

class Query {
  public:
    template <typename... Args>
    Query(const DeBruijnGraph &graph, Args&&... args)
        : graph_(graph), fwd_(std::forward<Args>(args)...) {
        if (std::any_of(fwd_.begin(), fwd_.end(), [](unsigned char c) { return c >= 128; }))
            throw std::invalid_argument("Alignment queries must contain ASCII characters");
#if ! _DNA_CASE_SENSITIVE_GRAPH
        // Match the column aligner's query normalization before creating views.
        std::transform(fwd_.begin(), fwd_.end(), fwd_.begin(),
                       [](unsigned char c) { return std::toupper(c); });
#endif
        rc_ = fwd_;
        reverse_complement(rc_.begin(), rc_.end());
    }

    const std::string& get_query(bool orientation = false) const {
        return !orientation ? fwd_ : rc_;
    }

    const DeBruijnGraph& get_graph() const { return graph_; }

  private:
    const DeBruijnGraph &graph_;
    std::string fwd_;
    std::string rc_;

};

} // namespace mtg::graph::align
