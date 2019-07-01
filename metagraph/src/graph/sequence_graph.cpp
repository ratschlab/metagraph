#include "sequence_graph.hpp"

#include <cassert>

typedef DeBruijnGraph::node_index node_index;


node_index DeBruijnGraph::kmer_to_node(const char *begin) const {
    return kmer_to_node(std::string(begin, get_k()));
}

node_index DeBruijnGraph::kmer_to_node(const std::string &kmer) const {
    assert(kmer.size() == get_k());

    node_index node = npos;
    map_to_nodes_sequentially(kmer.begin(), kmer.end(),
                              [&node](node_index i) { node = i; });
    return node;
}

// Check whether graph contains fraction of nodes from the sequence
bool DeBruijnGraph::find(const std::string &sequence,
                         double discovery_fraction) const {
    if (sequence.length() < get_k())
        return false;

    const size_t num_kmers = sequence.length() - get_k() + 1;
    const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    map_to_nodes(sequence,
        [&](node_index node) {
            if (node) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; }
    );

    return num_kmers_missing <= max_kmers_missing;
}

bool DeBruijnGraph::operator==(const DeBruijnGraph &) const {
    throw std::runtime_error("Not implemented");
    return false;
}

void DeBruijnGraph::call_nodes(const std::function<void(node_index)> &callback,
                               const std::function<bool()> &stop_early) const {
    const auto nnodes = num_nodes();
    for (node_index i = 1; i <= nnodes && !stop_early(); ++i) {
        callback(i);
    }
}

void call_sequences_from(const DeBruijnGraph &graph,
                         node_index start,
                         const std::function<void(const std::string&)> &callback,
                         sdsl::bit_vector *visited,
                         sdsl::bit_vector *discovered,
                         bool call_unitigs = false,
                         uint64_t min_tip_size = 0) {
    assert((min_tip_size <= 1 || call_unitigs)
                && "tip pruning works only for unitig extraction");
    assert(visited);
    assert(discovered);
    assert(!(*visited)[start]);

    std::stack<std::tuple<node_index, std::string, char>> paths;
    std::vector<std::pair<node_index, char>> targets;

    (*discovered)[start] = true;
    auto cur_node_seq = graph.get_node_sequence(start);
    paths.emplace(start,
                  cur_node_seq.substr(0, graph.get_k() - 1),
                  cur_node_seq.back());

    // keep traversing until we have worked off all branches from the queue
    while (paths.size()) {
        auto [node, sequence, next_char] = std::move(paths.top());
        paths.pop();
        start = node;

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (!(*visited)[node]) {
            assert(node);
            assert((*discovered)[node]);

            sequence.push_back(next_char);
            assert(sequence.length() >= graph.get_k());
            (*visited)[node] = true;

            targets.clear();
            graph.call_outgoing_kmers(node,
                [&](const auto &next, char c) { targets.emplace_back(next, c); }
            );

            if (targets.empty())
                break;

            if (targets.size() == 1
                    && (!call_unitigs || graph.indegree(targets[0].first) == 1)) {
                std::tie(node, next_char) = targets[0];
                (*discovered)[node] = true;
                continue;
            }

            auto new_seq = sequence.substr(sequence.length() - graph.get_k() + 1);
            node_index next_node = DeBruijnGraph::npos;
            //  _____.___
            //      \.___
            for (const auto& [next, c] : targets) {
                if (next_node == DeBruijnGraph::npos && !call_unitigs && !(*visited)[next]) {
                    (*discovered)[next] = true;
                    next_node = next;
                    next_char = c;
                } else if (!(*discovered)[next]) {
                    (*discovered)[next] = true;
                    paths.emplace(next, new_seq, c);
                }
            }

            if (next_node == DeBruijnGraph::npos)
                break;

            node = next_node;
        }

        if (sequence.size() >= graph.get_k()
              && (!call_unitigs
                  // check if not short
                  || sequence.size() >= graph.get_k() + min_tip_size - 1
                  // check if not tip
                  || graph.indegree(start) + graph.outdegree(node) >= 2)) {
            callback(sequence);
        }
    }
}

void DeBruijnGraph
::call_sequences(const std::function<void(const std::string&)> &callback) const {
    sdsl::bit_vector discovered(num_nodes() + 1, false);
    sdsl::bit_vector visited(num_nodes() + 1, false);

    // first, process start nodes (without incoming edges)
    call_nodes([&](const auto &start) {
        if (!visited[start] && !indegree(start)) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered);
        }
    });

    // then forks
    call_nodes([&](const auto &start) {
        if (!visited[start] && outdegree(start) > 1) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered);
        }
    });

    // then the rest of the cycles
    call_nodes([&](const auto &start) {
        if (!visited[start]) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered);
        }
    });
}

void DeBruijnGraph
::call_unitigs(const std::function<void(const std::string&)> &callback,
               size_t min_tip_size) const {
    sdsl::bit_vector discovered(num_nodes() + 1, false);
    sdsl::bit_vector visited(num_nodes() + 1, false);
    std::stack<std::tuple<node_index, node_index, std::string, char>> paths;
    std::vector<std::pair<node_index, char>> targets;

    // first, process start nodes (without incoming edges)
    call_nodes([&](const auto &start) {
        if (!visited[start] && !indegree(start)) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered,
                                true,
                                min_tip_size);
        }
    });

    // then merges
    //  ____.____
    //  ___/
    //
    call_nodes([&](const auto &start) {
        if (!visited[start] && indegree(start) > 1) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered,
                                true,
                                min_tip_size);
        }
    });

    // then forks
    //  ____.____
    //       \___
    //
    call_nodes([&](const auto &start) {
        if (!visited[start] && outdegree(start) > 1) {
            call_outgoing_kmers(start, [&](node_index next, char) {
                if (!visited[next])
                    call_sequences_from(*this,
                                        next,
                                        callback,
                                        &visited,
                                        &discovered,
                                        true,
                                        min_tip_size);
            });
        }
    });

    // then the rest of the cycles
    call_nodes([&](const auto &start) {
        if (!visited[start]) {
            call_sequences_from(*this,
                                start,
                                callback,
                                &visited,
                                &discovered,
                                true,
                                min_tip_size);
        }
    });
}

/**
 * Traverse graph and iterate over all nodes
 */
void DeBruijnGraph
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    sdsl::bit_vector visited(num_nodes() + 1, false);
    std::stack<std::pair<node_index, std::string>> nodes;

    call_nodes([&](node_index i) {
        if (visited[i])
            return;

        nodes.emplace(i, get_node_sequence(i));
        while (nodes.size()) {
            // FYI: structured binding is a new thing that often
            // leads to issues, so avoid using it.
            // https://stackoverflow.com/questions/50799719/reference-to-local-binding-declared-in-enclosing-function?noredirect=1&lq=1
            auto [node, sequence] = std::move(nodes.top());
            nodes.pop();

            while (!visited[node]) {
                visited[node] = true;
                callback(node, sequence);
                sequence.assign(sequence.begin() + 1, sequence.end());

                auto next_node = npos;
                char next_c = '\0';
                call_outgoing_kmers(
                    node,
                    [&, &sequence=sequence](const auto &next, char c) {
                        if (visited[next])
                            return;

                        if (next_node == npos) {
                            std::tie(next_node, next_c) =  std::tie(next, c);
                        } else {
                            nodes.emplace(next, sequence + c);
                        }
                    }
                );

                if (next_node == npos)
                    break;

                node = next_node;
                sequence.push_back(next_c);
            }
        }
    });
}

void DeBruijnGraph::print(std::ostream &out) const {
    auto vertex_header = std::string("Vertex");
    vertex_header.resize(get_k(), ' ');

    out << "Index"
        << "\t" << vertex_header
        << std::endl;

    call_nodes([&](node_index i) {
        out << i
            << "\t" << get_node_sequence(i)
            << std::endl;
    });
}

std::ostream& operator<<(std::ostream &out, const DeBruijnGraph &graph) {
    graph.print(out);
    return out;
}
