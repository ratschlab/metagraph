#include "sequence_graph.hpp"

#include <cassert>
#include <progress_bar.hpp>

#include "threading.hpp"
#include "bitmap.hpp"

namespace utils {
    bool get_verbose();
}

typedef DeBruijnGraph::node_index node_index;

static const uint64_t kBlockSize = 9'999'872;
static_assert(!(kBlockSize & 0xFF));


/*************** SequenceGraph ***************/

void SequenceGraph::add_extension(std::shared_ptr<GraphExtension> extension) {
    assert(extension.get());
    extensions_.push_back(extension);
}

void SequenceGraph::serialize_extensions(const std::string &filename) const {
    for (auto extension : extensions_) {
        extension->serialize(utils::remove_suffix(filename, file_extension()) + file_extension());
    }
}

/*************** DeBruijnGraph ***************/

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

void DeBruijnGraph::traverse(node_index start,
                             const char* begin,
                             const char* end,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    assert(in_graph(start));
    assert(end >= begin);
    if (terminate())
        return;

    for (; begin != end && !terminate(); ++begin) {
        start = traverse(start, *begin);

        if (start == npos)
            return;

        callback(start);
    }
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
                         ProgressBar &progress_bar,
                         bool call_unitigs = false,
                         uint64_t min_tip_size = 0) {
    assert(graph.in_graph(start));
    assert((min_tip_size <= 1 || call_unitigs)
                && "tip pruning works only for unitig extraction");
    assert(visited);
    assert(discovered);
    assert(!(*visited)[start]);

    std::vector<node_index> queue = { start };
    (*discovered)[start] = true;

    std::vector<std::pair<node_index, char>> targets;

    // keep traversing until we have worked off all branches from the queue
    while (queue.size()) {
        node_index node = queue.back();
        std::string sequence = graph.get_node_sequence(node);
        queue.pop_back();
        start = node;
        assert(!call_unitigs || !(*visited)[node]);
        if ((*visited)[node])
            continue;

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (!(*visited)[node]) {
            assert(node);
            assert((*discovered)[node]);
            assert(sequence.length() >= graph.get_k());
            (*visited)[node] = true;
            ++progress_bar;

            targets.clear();
            graph.call_outgoing_kmers(node,
                [&](auto next, char c) { targets.emplace_back(next, c); }
            );

            if (targets.empty())
                break;

            // in call_unitigs mode, all nodes with multiple incoming
            // edges are marked as discovered
            assert(!call_unitigs || graph.has_single_incoming(targets.front().first)
                                 || (*discovered)[targets.front().first]);
            if (targets.size() == 1) {
                if ((*visited)[targets.front().first])
                    break;

                if (!call_unitigs || !(*discovered)[targets.front().first]) {
                    sequence.push_back('\0');
                    std::tie(node, sequence.back()) = targets[0];
                    (*discovered)[node] = true;
                    continue;
                }
            }

            node_index next_node = DeBruijnGraph::npos;
            //  _____.___
            //      \.___
            for (const auto& [next, c] : targets) {
                if (next_node == DeBruijnGraph::npos
                        && !call_unitigs
                        && !(*visited)[next]) {
                    (*discovered)[next] = true;
                    next_node = next;
                    sequence.push_back(c);
                } else if (!(*discovered)[next]) {
                    (*discovered)[next] = true;
                    queue.push_back(next);
                }
            }

            if (next_node == DeBruijnGraph::npos)
                break;

            node = next_node;
        }

        assert(sequence.size() >= graph.get_k());

        if (!call_unitigs
                  // check if long
                  || sequence.size() >= graph.get_k() + min_tip_size - 1
                  // check if not tip
                  || graph.indegree(start) + graph.outdegree(node) >= 2) {
            callback(sequence);
        }
    }
}

void call_sequences(const DeBruijnGraph &graph,
                    const std::function<void(const std::string&)> &callback,
                    bool call_unitigs,
                    uint64_t min_tip_size = 0) {
    sdsl::bit_vector discovered(graph.num_nodes() + 1, true);
    graph.call_nodes([&](auto node) { discovered[node] = false; });
    sdsl::bit_vector visited = discovered;

    ProgressBar progress_bar(discovered.size() - sdsl::util::cnt_one_bits(visited),
                             "Traverse graph",
                             std::cerr, !utils::get_verbose());

    auto call_paths_from = [&](node_index node) {
        call_sequences_from(graph,
                            node,
                            callback,
                            &visited,
                            &discovered,
                            progress_bar,
                            call_unitigs,
                            min_tip_size);
    };

    if (call_unitigs) {
        // mark all source and merge nodes (those with indegree 0 or >1)
        //  .____  or  .____  or  ____.___
        //              \___      ___/
        //
        #pragma omp parallel for num_threads(get_num_threads())
        for (uint64_t begin = 0; begin < discovered.size(); begin += kBlockSize) {
            call_zeros(discovered,
                begin,
                std::min(begin + kBlockSize, discovered.size()),
                [&](auto i) { discovered[i] = !graph.has_single_incoming(i); }
            );
        }

        // now traverse graph starting at these nodes
        call_zeros(visited, [&](auto node) {
            assert(discovered[node] == !graph.has_single_incoming(node));
            if (discovered[node])
                call_paths_from(node);
        });

    } else {
        // start at the source nodes (those with indegree == 0)
        //  .____  or  .____
        //              \___
        //
        graph.call_source_nodes([&](auto node) {
            assert(!visited[node]);
            assert(graph.has_no_incoming(node));

            call_paths_from(node);
        });
    }

    // then forks
    //  ____.____
    //       \___
    //
    call_zeros(visited, [&](auto node) {
        if (graph.has_multiple_outgoing(node)) {
            graph.adjacent_outgoing_nodes(node, [&](auto next) {
                if (!visited[next])
                    call_paths_from(next);
            });
        }
    });

    // then the rest (loops)
    call_zeros(visited, call_paths_from);
}

void DeBruijnGraph
::call_sequences(const std::function<void(const std::string&)> &callback) const {
    ::call_sequences(*this, callback, false);
}

void DeBruijnGraph
::call_unitigs(const std::function<void(const std::string&)> &callback,
               size_t min_tip_size) const {
    ::call_sequences(*this, callback, true, min_tip_size);
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
            // https://bit.ly/2JI6oci
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

void DeBruijnGraph
::call_source_nodes(const std::function<void(node_index)> &callback) const {
    call_nodes([&](node_index i) {
        if (has_no_incoming(i))
            callback(i);
    });
}


std::ostream& operator<<(std::ostream &out, const DeBruijnGraph &graph) {
    graph.print(out);
    return out;
}

// returns the edge rank, starting from zero
size_t incoming_edge_rank(const DeBruijnGraph &graph,
                          DeBruijnGraph::node_index source,
                          DeBruijnGraph::node_index target) {
    assert(graph.in_graph(source));
    assert(graph.in_graph(target));

    assert(graph.get_node_sequence(source).substr(1)
                == graph.get_node_sequence(target).substr(0, graph.get_k() - 1));

    size_t edge_rank = 0;
    bool done = false;

    graph.adjacent_incoming_nodes(target, [&](auto node) {
        if (node == source)
            done = true;

        if (!done)
            edge_rank++;
    });

    assert(done && "the edge must exist in graph");

    return edge_rank;
}
