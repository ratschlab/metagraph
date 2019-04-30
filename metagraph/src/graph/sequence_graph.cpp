#include "sequence_graph.hpp"

#include <cassert>


DeBruijnGraph::node_index DeBruijnGraph::kmer_to_node(const char *begin) const {
    return kmer_to_node(std::string(begin, get_k()));
}

DeBruijnGraph::node_index DeBruijnGraph::kmer_to_node(const std::string &kmer) const {
    assert(kmer.size() == get_k());

    node_index node = npos;
    map_to_nodes(kmer, [&node](node_index i) { node = i; });
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

void DeBruijnGraph::call_nodes(const std::function<void(const node_index&)> &callback,
                               const std::function<bool()> &stop_early) const {
    auto nnodes = num_nodes();
    for (node_index i = 1; i <= nnodes; ++i) {
        callback(i);
        if (stop_early())
            return;
    }
}

void DeBruijnGraph
::call_sequences(const std::function<void(const std::string&)> &callback) const {
    auto nnodes = num_nodes();
    sdsl::bit_vector discovered(nnodes + 1, false);
    sdsl::bit_vector visited(nnodes + 1, false);
    std::stack<std::tuple<node_index, node_index, std::string, char>> paths;
    std::vector<std::pair<node_index, char>> targets;

    call_nodes(
        [&](const auto &start) {
            if (visited[start])
                return;

            call_sequences_from(start,
                                callback,
                                &visited,
                                &discovered,
                                &paths,
                                &targets);
            }
    );
}

void DeBruijnGraph
::call_unitigs(const std::function<void(const std::string&)> &callback,
               size_t max_pruned_dead_end_size) const {
    auto nnodes = num_nodes();
    sdsl::bit_vector discovered(nnodes + 1, false);
    sdsl::bit_vector visited(nnodes + 1, false);
    std::stack<std::tuple<node_index, node_index, std::string, char>> paths;
    std::vector<std::pair<node_index, char>> targets;

    call_nodes(
        [&](const auto &start) {
            if (visited[start])
                return;

            call_sequences_from(start,
                                callback,
                                &visited,
                                &discovered,
                                &paths,
                                &targets,
                                true,
                                max_pruned_dead_end_size);
            }
    );
}


void DeBruijnGraph
::call_sequences_from(node_index start,
                      const std::function<void(const std::string&)> &callback,
                      sdsl::bit_vector *visited,
                      sdsl::bit_vector *discovered,
                      std::stack<std::tuple<node_index, node_index, std::string, char>> *paths,
                      std::vector<std::pair<node_index, char>> *targets,
                      bool split_to_contigs,
                      uint64_t max_pruned_dead_end_size) const {
    assert(visited);
    assert(discovered);
    assert(paths);
    assert(targets);
    assert(!(*visited)[start]);

    (*discovered)[start] = true;
    auto cur_node_seq = get_node_sequence(start);
    paths->emplace(start,
                   start,
                   cur_node_seq.substr(0, get_k() - 1),
                   cur_node_seq.back());

    // keep traversing until we have worked off all branches from the queue
    while (paths->size()) {
        auto [first, node, sequence, next_char] = std::move(paths->top());
        paths->pop();

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (!(*visited)[node]) {
            assert(node);
            assert((*discovered)[node]);

            sequence.push_back(next_char);
            assert(sequence.length() >= get_k());
            (*visited)[node] = true;

            targets->clear();
            call_outgoing_kmers(node,
                                [&](const auto &next, char c) {
                                    targets->emplace_back(next, c);
                                });

            if (targets->empty())
                break;

            bool continue_traversal = !split_to_contigs || indegree(node) == 1;

            if (continue_traversal && targets->size() == 1) {
                std::tie(node, next_char) = (*targets)[0];
                (*discovered)[node] = true;
                continue;
            }

            auto new_seq = sequence.substr(sequence.length() - get_k() + 1);
            node_index next_node = npos;
            for (const auto& [next, c] : *targets) {
                if (next_node == npos && !split_to_contigs && !(*visited)[next]) {
                    (*discovered)[next] = true;
                    next_node = next;
                    next_char = c;
                } else if (!(*discovered)[next]) {
                    (*discovered)[next] = true;
                    paths->emplace(node, next, new_seq, c);
                }
            }

            if (next_node == npos)
                break;

            node = next_node;
        }
        if (sequence.size() >= get_k()
              && (!split_to_contigs
              || indegree(first)
              || outdegree(node)
              || sequence.size() >= get_k() + max_pruned_dead_end_size))
            callback(sequence);
    }
}

/**
 * Traverse graph and iterate over all nodes
 */
void DeBruijnGraph
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    auto nnodes = num_nodes();
    sdsl::bit_vector visited(nnodes + 1, false);
    std::stack<std::pair<node_index, std::string>> nodes;

    for (node_index i = 1; i <= nnodes; ++i) {
        if (visited[i])
            continue;

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
    }
}
