#include "aln_match.hpp"

#include "common/logger.hpp"

namespace mtg::graph::align_redone {

using mtg::common::logger;

// Return the string spelled by the path. This path may have disconnects (if it came)
// from a chain alignment), so this method handles that case. If there is an invalid
// edge, or if there is too long of a stretch of npos nodes, this throws a runtime error.
std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path) {
    std::string seq;

    if (path.empty())
        return seq;

    seq.reserve(path.size() + graph.get_k() - 1);

    size_t num_dummy = 0;
    size_t num_unknown = 0;
    if (path.front()) {
        seq += graph.get_node_sequence(path.front());
    } else {
        num_unknown = graph.get_k();
        seq += std::string(num_unknown, '$');
        num_dummy = 1;
    }

    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i]) {
            if (num_dummy) {
                seq += '$';
                ++num_unknown;
                std::string next_seq = graph.get_node_sequence(path[i]);
                auto it = seq.end() - next_seq.size();
                for (char c : next_seq) {
                    if (*it == '$' && c != '$') {
                        --num_unknown;
                        *it = c;
                    }

                    ++it;
                }
                num_dummy = 0;
            } else {
                char next = '\0';
                graph.call_outgoing_kmers(path[i - 1], [&](auto next_node, char c) {
                    if (next_node == path[i])
                        next = c;
                });

                if (!next) {
                    logger->error("Invalid edge {} {}\t{} {}",
                                  path[i - 1], path[i],
                                  graph.get_node_sequence(path[i - 1]),
                                  graph.get_node_sequence(path[i]));
                    throw std::runtime_error("");
                }

                seq += next;
                if (num_unknown) {
                    std::string next_seq = graph.get_node_sequence(path[i]);
                    auto it = seq.end() - next_seq.size();
                    for (char c : next_seq) {
                        if (*it == '$' && c != '$') {
                            --num_unknown;
                            *it = c;
                        }

                        ++it;
                    }
                }
            }
        } else {
            seq += '$';
            ++num_dummy;
            ++num_unknown;
        }
    }

    assert(seq.size() == path.size() + graph.get_k() - 1);
    return seq;
}

void Anchor::append(const Anchor &other, const DeBruijnGraph *graph) {
    assert(query_.data() == other.query_.data());
    assert(query_.data() + query_.size() == other.query_.data() + other.query_.size());
    assert(label_class_ == other.label_class_);

    if (empty()) {
        *this = other;
        return;
    }

    assert(other.seed_.data() >= seed_.data());

    if (seed_.data() + seed_.size() >= other.seed_.data() + other.seed_.size())
        return;

    if (seed_.data() == other.seed_.data()) {
        assert(coord_ == other.coord_);
        *this = other;
        return;
    }

    auto next = path_.begin() + (other.seed_.data() - seed_.data());
    assert(next <= path_.end());

    seed_ = std::string_view(
        seed_.data(),
        other.seed_.data() + other.seed_.size() - seed_.data()
    );

    path_.erase(next, path_.end());
    path_.insert(path_.end(), other.path_.begin(), other.path_.end());
    suffix_ = other.suffix_;

    assert(!graph || is_spelling_valid(*graph));
}

} // namespace mtg::graph::align