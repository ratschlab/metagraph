#pragma once

#include <string>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg::graph::align {

class Query {
  public:
    template <typename... Args>
    Query(const DeBruijnGraph &graph, Args&&... args)
        : graph_(graph), fwd_(std::forward<Args>(args)...), rc_(generate_rc()) {}

    const std::string& get_query(bool orientation = false) const {
        return !orientation ? fwd_ : rc_;
    }

    const DeBruijnGraph& get_graph() const { return graph_; }

  private:
    const DeBruijnGraph &graph_;
    std::string fwd_;
    std::string rc_;

    std::string generate_rc() {
        std::string rc = fwd_;
        reverse_complement(rc_.begin(), rc_.end());
        return rc;
    }
};

} // namespace mtg::graph::align