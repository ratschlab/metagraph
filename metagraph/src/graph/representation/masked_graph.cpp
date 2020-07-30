#include "masked_graph.hpp"

#include "common/serialization.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                      std::unique_ptr<bitmap>&& kmers_in_graph,
                      bool only_valid_nodes_in_mask,
                      bool canonical)
      : graph_(graph),
        kmers_in_graph_(std::move(kmers_in_graph)),
        only_valid_nodes_in_mask_(only_valid_nodes_in_mask),
        is_canonical_(canonical) {
    assert(kmers_in_graph_.get());
    assert(kmers_in_graph_->size() == graph->max_index() + 1);
    assert(!is_canonical_ || graph_->is_canonical_mode());
}

MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                      std::function<bool(node_index)>&& callback,
                      bool only_valid_nodes_in_mask,
                      bool canonical)
      : MaskedDeBruijnGraph(graph,
                            std::make_unique<bitmap_lazy>(std::move(callback),
                                                          graph->max_index() + 1),
                            only_valid_nodes_in_mask,
                            canonical) {}

// Traverse the outgoing edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse(node_index node, char next_char) const {
    assert(in_subgraph(node));

    auto index = graph_->traverse(node, next_char);
    return index && in_subgraph(index) ? index : DeBruijnGraph::npos;
}
// Traverse the incoming edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse_back(node_index node, char prev_char) const {
    assert(in_subgraph(node));

    auto index = graph_->traverse_back(node, prev_char);
    return index && in_subgraph(index) ? index : DeBruijnGraph::npos;
}

size_t MaskedDeBruijnGraph::outdegree(node_index node) const {
    assert(in_subgraph(node));

    size_t outdegree = 0;
    graph_->adjacent_outgoing_nodes(node, [&](auto index) { outdegree += in_subgraph(index); });
    return outdegree;
}

size_t MaskedDeBruijnGraph::indegree(node_index node) const {
    assert(in_subgraph(node));

    size_t indegree = 0;
    graph_->adjacent_incoming_nodes(node, [&](auto index) { indegree += in_subgraph(index); });
    return indegree;
}

void MaskedDeBruijnGraph
::adjacent_outgoing_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_subgraph(node));

    graph_->adjacent_outgoing_nodes(node, [&](auto node) {
        if (in_subgraph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::adjacent_incoming_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_subgraph(node));

    graph_->adjacent_incoming_nodes(node, [&](auto node) {
        if (in_subgraph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::call_outgoing_kmers(node_index kmer,
                      const OutgoingEdgeCallback &callback) const {
    assert(in_subgraph(kmer));

    graph_->call_outgoing_kmers(
        kmer,
        [&](const auto &index, auto c) {
            if (in_subgraph(index))
                callback(index, c);
        }
    );
}

void MaskedDeBruijnGraph
::call_incoming_kmers(node_index kmer,
                      const IncomingEdgeCallback &callback) const {
    assert(in_subgraph(kmer));

    graph_->call_incoming_kmers(
        kmer,
        [&](const auto &index, auto c) {
            if (in_subgraph(index))
                callback(index, c);
        }
    );
}

bit_vector_stat get_boss_mask(const DBGSuccinct &dbg_succ,
                              const bitmap &kmers_in_graph,
                              bool only_valid_nodes_in_mask) {
    sdsl::bit_vector mask_bv(dbg_succ.get_boss().num_edges() + 1, false);
    if (only_valid_nodes_in_mask) {
        kmers_in_graph.call_ones(
            [&](auto i) {
                assert(dbg_succ.kmer_to_boss_index(i));
                mask_bv[dbg_succ.kmer_to_boss_index(i)] = true;
            }
        );
    } else {
        dbg_succ.call_nodes(
            [&](auto i) {
                assert(dbg_succ.kmer_to_boss_index(i));
                if (kmers_in_graph[i])
                    mask_bv[dbg_succ.kmer_to_boss_index(i)] = true;
            }
        );
    }
    return bit_vector_stat(std::move(mask_bv));
}

void MaskedDeBruijnGraph
::call_sequences(const CallPath &callback,
                 size_t num_threads,
                 bool kmers_in_single_form) const {
    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {
        bit_vector_stat mask = get_boss_mask(*dbg_succ, *kmers_in_graph_,
                                             only_valid_nodes_in_mask_);

        dbg_succ->get_boss().call_sequences([&](std::string&& sequence, auto&& path) {
            for (auto &node : path) {
                node = dbg_succ->boss_to_kmer_index(node);
            }
            callback(sequence, path);

        }, num_threads, kmers_in_single_form, &mask);

    } else {
        DeBruijnGraph::call_sequences(callback, num_threads, kmers_in_single_form);
    }
}

void MaskedDeBruijnGraph
::call_unitigs(const CallPath &callback,
               size_t num_threads,
               size_t min_tip_size,
               bool kmers_in_single_form) const {
    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {
        bit_vector_stat mask = get_boss_mask(*dbg_succ, *kmers_in_graph_,
                                             only_valid_nodes_in_mask_);

        dbg_succ->get_boss().call_unitigs([&](std::string&& sequence, auto&& path) {
            for (auto &node : path) {
                node = dbg_succ->boss_to_kmer_index(node);
            }
            callback(sequence, path);

        }, num_threads, min_tip_size, kmers_in_single_form, &mask);

    } else {
        DeBruijnGraph::call_unitigs(callback,
                                    num_threads,
                                    min_tip_size,
                                    kmers_in_single_form);
    }
}

void MaskedDeBruijnGraph
::call_nodes(const std::function<void(node_index)> &callback,
             const std::function<bool()> &stop_early) const {
    assert(max_index() + 1 == kmers_in_graph_->size());

    bool stop = false;

    if (only_valid_nodes_in_mask_) {
        // iterate only through the nodes marked in the mask
        // TODO: add terminate<bool(void)> to call_ones
        kmers_in_graph_->call_ones(
            [&](auto index) {
                if (stop || !index)
                    return;

                assert(in_subgraph(index));

                if (stop_early()) {
                    stop = true;
                } else {
                    callback(index);
                }
            }
        );
    } else {
        // call all nodes in the base graph and check the mask
        graph_->call_nodes(
            [&](auto index) {
                if (in_subgraph(index))
                    callback(index);
            },
            stop_early
        );
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
void MaskedDeBruijnGraph
::map_to_nodes(std::string_view sequence,
               const std::function<void(node_index)> &callback,
               const std::function<bool()> &terminate) const {
    graph_->map_to_nodes(
        sequence,
        [&](const node_index &index) {
            callback(index && in_subgraph(index) ? index : npos);
        },
        terminate
    );
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence
void MaskedDeBruijnGraph
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    graph_->map_to_nodes_sequentially(
        sequence,
        [&](const node_index &index) {
            callback(index && in_subgraph(index) ? index : npos);
        },
        terminate
    );
}

// Get string corresponding to |node_index|.
// Note: Not efficient if sequences in nodes overlap. Use sparingly.
std::string MaskedDeBruijnGraph::get_node_sequence(node_index index) const {
    assert(in_subgraph(index));

    return graph_->get_node_sequence(index);
}

bool MaskedDeBruijnGraph::operator==(const MaskedDeBruijnGraph &other) const {
    return get_k() == other.get_k()
            && is_canonical_mode() == other.is_canonical_mode()
            && num_nodes() == other.num_nodes()
            && *kmers_in_graph_ == *other.kmers_in_graph_
            && *graph_ == *other.graph_;
}

bool MaskedDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (dynamic_cast<const MaskedDeBruijnGraph*>(&other))
        return operator==(dynamic_cast<const MaskedDeBruijnGraph&>(other));

    return DeBruijnGraph::operator==(other);
}

void MaskedDeBruijnGraph::update_mask(const GenerateNodes &generate_nodes,
                                      bool only_valid_nodes_in_mask,
                                      bool canonical,
                                      bool async,
                                      int memorder) {
    only_valid_nodes_in_mask_ = only_valid_nodes_in_mask;
    is_canonical_ = canonical;

    auto *mask_updateable = dynamic_cast<bitmap_dyn*>(kmers_in_graph_.get());
    if (mask_updateable) {
        auto *mask_vector = dynamic_cast<bitmap_vector*>(mask_updateable);
        if (mask_vector) {
            auto &mask_sdsl = const_cast<sdsl::bit_vector&>(mask_vector->data());
            generate_nodes([&](node_index i, bool val) {
                if (val) {
                    set_bit(mask_sdsl.data(), i, async, memorder);
                } else {
                    unset_bit(mask_sdsl.data(), i, async, memorder);
                }
            });
        } else {
            generate_nodes([&](node_index i, bool val) { mask_updateable->set(i, val); });
        }
    } else {
        sdsl::bit_vector mask_updated(kmers_in_graph_->size(), false);
        kmers_in_graph_->add_to(&mask_updated);
        generate_nodes([&](node_index i, bool val) { mask_updated[i] = val; });

        // if kmers_in_graph_ was not bitmap_dyn before, don't make it now
        kmers_in_graph_ = std::make_unique<bit_vector_stat>(std::move(mask_updated));
    }
}

} // namespace graph
} // namespace mtg
