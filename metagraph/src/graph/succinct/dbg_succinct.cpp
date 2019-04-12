#include "dbg_succinct.hpp"

#include <cassert>
#include <vector>
#include <algorithm>
#include <string>

#include "boss.hpp"
#include "serialization.hpp"
#include "reverse_complement.hpp"
#include "utils.hpp"

using utils::remove_suffix;

typedef DBGSuccinct::node_index node_index;


DBGSuccinct::DBGSuccinct(size_t k, bool canonical_mode)
      : boss_graph_(std::make_unique<BOSS>(k - 1)),
        canonical_mode_(canonical_mode) {}

DBGSuccinct::DBGSuccinct(BOSS *boss_graph, bool canonical_mode)
      : boss_graph_(boss_graph),
        canonical_mode_(canonical_mode) {}

size_t DBGSuccinct::get_k() const {
    return boss_graph_->get_k() + 1;
}

// Check whether graph contains fraction of nodes from the sequence
bool DBGSuccinct::find(const std::string &sequence,
                       double discovery_fraction) const {
    return boss_graph_->find(sequence, discovery_fraction);
}

// Traverse the outgoing edge
node_index DBGSuccinct::traverse(node_index node, char next_char) const {
    assert(node);

    // dbg node is a boss edge
    BOSS::edge_index edge = boss_graph_->fwd(kmer_to_boss_index(node));
    return boss_to_kmer_index(
        boss_graph_->pick_edge(edge,
                               boss_graph_->get_source_node(edge),
                               boss_graph_->encode(next_char))
    );
}

// Traverse the incoming edge
node_index DBGSuccinct::traverse_back(node_index node, char prev_char) const {
    assert(node);

    // map dbg node, i.e. a boss edge, to a boss node
    auto boss_edge = kmer_to_boss_index(node);

    auto boss_node = boss_graph_->get_source_node(boss_edge);
    auto source_node = boss_graph_->traverse_back(boss_node, prev_char);

    if (!source_node)
        return npos;

    return boss_to_kmer_index(
        boss_graph_->outgoing_edge_idx(source_node,
                                       boss_graph_->get_node_last_value(boss_edge))
    );
}

void DBGSuccinct::call_outgoing_kmers(node_index node,
                                      const OutgoingEdgeCallback &callback) const {
    assert(node);

    auto boss_edge = kmer_to_boss_index(node);

    auto last = boss_graph_->fwd(boss_edge);
    auto first = boss_graph_->pred_last(last - 1) + 1;

    for (auto i = first; i <= last; ++i) {
        assert(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size
                == boss_graph_->get_node_last_value(i));

        auto next = boss_to_kmer_index(i);
        if (next != npos)
            callback(next, boss_graph_->decode(boss_graph_->get_W(i)));
    }
}

void DBGSuccinct::adjacent_outgoing_nodes(node_index node,
                                          std::vector<node_index> *target_nodes) const {
    assert(node);
    assert(target_nodes);

    auto boss_edge = kmer_to_boss_index(node);

    auto last = boss_graph_->fwd(boss_edge);
    auto first = boss_graph_->pred_last(last - 1) + 1;

    for (auto i = first; i <= last; ++i) {
        assert(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size
                == boss_graph_->get_node_last_value(i));

        auto next = boss_to_kmer_index(i);
        if (next != npos)
            target_nodes->emplace_back(next);
    }
}

void DBGSuccinct::adjacent_incoming_nodes(node_index node,
                                          std::vector<node_index> *source_nodes) const {
    assert(node);
    assert(source_nodes);

    auto edge = kmer_to_boss_index(node);

    boss_graph_->call_adjacent_incoming_edges(edge,
        [&](BOSS::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = boss_to_kmer_index(incoming_boss_edge);
            if (prev != npos)
                source_nodes->emplace_back(prev);
        }
    );
}

// Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
// is passed. If passed, |nodes_inserted| must have length equal
// to the number of nodes in graph.
void DBGSuccinct::add_sequence(const std::string &sequence,
                               bit_vector_dyn *nodes_inserted) {
    add_seq(sequence, nodes_inserted);
    if (canonical_mode_) {
        // insert reverse complement sequence as well,
        // to have all canonical k-mers in graph
        std::string sequence_copy = sequence;
        reverse_complement(sequence_copy.begin(), sequence_copy.end());
        add_seq(sequence_copy, nodes_inserted);
    }
}

void DBGSuccinct::add_seq(const std::string &sequence,
                          bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    if (nodes_inserted) {
        std::vector<uint64_t> inserted_indexes;
        inserted_indexes.reserve(sequence.size());

        boss_graph_->add_sequence(sequence, true, &inserted_indexes);

        for (auto i : inserted_indexes) {
            if (valid_edges_.get())
                valid_edges_->insert_bit(i, true);
            nodes_inserted->insert_bit(boss_to_kmer_index(i), true);
        }
    } else {
        boss_graph_->add_sequence(sequence, true);
    }

    assert(!valid_edges_.get() || !(*valid_edges_)[0]);
}

std::string DBGSuccinct::get_node_sequence(node_index node) const {
    assert(node);
    assert(node <= num_nodes());

    auto boss_edge = kmer_to_boss_index(node);

    return boss_graph_->get_node_str(boss_edge)
            + boss_graph_->decode(boss_graph_->get_W(boss_edge));
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are not mapped to canonical ones
void DBGSuccinct::map_to_nodes_sequentially(std::string::const_iterator begin,
                                            std::string::const_iterator end,
                                            const std::function<void(node_index)> &callback,
                                            const std::function<bool()> &terminate) const {
    if (begin + get_k() > end)
        return;

    boss_graph_->map_to_edges(
        std::string(begin, end),
        [&](BOSS::edge_index i) { callback(boss_to_kmer_index(i)); },
        terminate
    );
}

// Map sequence k-mers to the canonical graph nodes
// and run callback for each node until the termination condition is satisfied
void DBGSuccinct::map_to_nodes(const std::string &sequence,
                               const std::function<void(node_index)> &callback,
                               const std::function<bool()> &terminate) const {
    if (sequence.size() < get_k())
        return;

    if (canonical_mode_) {
        auto forward = boss_graph_->map_to_edges(sequence);

        std::string sequence_rev_compl = sequence;
        reverse_complement(sequence_rev_compl.begin(), sequence_rev_compl.end());

        auto rev_compl = boss_graph_->map_to_edges(sequence_rev_compl);

        assert(forward.size() == sequence.size() - get_k() + 1);
        assert(forward.size() == rev_compl.size());

        for (size_t i = 0; i < forward.size() && !terminate(); ++i) {
            // the definition of a canonical k-mer is redefined:
            //      use k-mer with smaller index in the BOSS table.
            callback(boss_to_kmer_index(
                std::min(forward[i], rev_compl[rev_compl.size() - 1 - i])
            ));
        }

    } else {
        boss_graph_->map_to_edges(
            sequence,
            [&](BOSS::edge_index i) { callback(boss_to_kmer_index(i)); },
            terminate
        );
    }
}

size_t DBGSuccinct::outdegree(node_index node) const {
    assert(node);

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge > 1 && !boss_graph_->get_W(boss_edge)) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return 0;
    }

    auto last_target_kmer = boss_graph_->fwd(boss_edge);

    if (!boss_graph_->get_W(last_target_kmer)) {
        // There is a sunk dummy target, hence this is the only outgoing edge
        if (valid_edges_.get()) {
            // if mask is defined, skip boss dummy sink edges
            return 0;
        } else {
            return 1;
        }
    }

    return last_target_kmer - boss_graph_->pred_last(last_target_kmer - 1);
}

size_t DBGSuccinct::indegree(node_index node) const {
    assert(node);

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return 1;

    auto d = boss_graph_->get_node_last_value(boss_edge);

    auto x = boss_graph_->bwd(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    if (x + 1 == boss_graph_->get_W().size())
        return first_valid;

    uint64_t y = boss_graph_->succ_W(x + 1, d);
    return first_valid + boss_graph_->rank_W(y - 1, d + boss_graph_->alph_size)
                        - boss_graph_->rank_W(x - 1, d + boss_graph_->alph_size);
}

uint64_t DBGSuccinct::num_nodes() const {
    return valid_edges_.get()
                ? valid_edges_->num_set_bits()
                : boss_graph_->num_edges();
}

bool DBGSuccinct::load(const std::string &filename) {
    {
        std::ifstream instream(remove_suffix(filename, kExtension) + kExtension,
                               std::ios::binary);

        if (!boss_graph_->load(instream))
            return false;

        try {
            canonical_mode_ = load_number(instream);
        } catch (...) {
            canonical_mode_ = false;
        }
    }

    // release the old mask
    valid_edges_.reset();

    std::ifstream instream(remove_suffix(filename, kExtension) + kDummyMaskExtension,
                           std::ios::binary);
    if (!instream.good())
        return true;

    // initialize a new vector
    switch (get_state()) {
        case Config::STAT: {
            valid_edges_.reset(new bit_vector_stat());
            break;
        }
        case Config::DYN: {
            valid_edges_.reset(new bit_vector_dyn());
            break;
        }
        case Config::SMALL: {
            valid_edges_.reset(new bit_vector_small());
            break;
        }
    }

    // load the mask of valid edges (all non-dummy including npos 0)
    if (!valid_edges_->load(instream)) {
        std::cerr << "Error: Can't load dummy edge mask." << std::endl;
        return false;
    }

    if (valid_edges_->size() != boss_graph_->num_edges() + 1 || (*valid_edges_)[0]) {
        std::cerr << "Error: Edge mask is not compatible with graph." << std::endl;
        return false;
    }

    return true;
}

void DBGSuccinct::serialize(const std::string &filename) const {
    {
        const auto out_filename = remove_suffix(filename, kExtension) + kExtension;
        std::ofstream outstream(out_filename, std::ios::binary);
        boss_graph_->serialize(outstream);
        serialize_number(outstream, canonical_mode_);

        if (!outstream.good())
            throw std::ios_base::failure("Can't write to file " + out_filename);
    }

    if (!valid_edges_.get())
        return;

    assert((boss_graph_->get_state() == Config::StateType::STAT
                && dynamic_cast<const bit_vector_stat*>(valid_edges_.get()))
        || (boss_graph_->get_state() == Config::StateType::DYN
                && dynamic_cast<const bit_vector_dyn*>(valid_edges_.get()))
        || (boss_graph_->get_state() == Config::StateType::SMALL
                && dynamic_cast<const bit_vector_small*>(valid_edges_.get())));

    const auto out_filename = remove_suffix(filename, kExtension) + kDummyMaskExtension;
    std::ofstream outstream(out_filename, std::ios::binary);
    if (!outstream.good())
        throw std::ios_base::failure("Can't write to file " + out_filename);

    valid_edges_->serialize(outstream);
}

void DBGSuccinct::switch_state(Config::StateType new_state) {
    if (get_state() == new_state)
        return;

    if (valid_edges_.get()) {
        switch (new_state) {
            case Config::STAT: {
                valid_edges_ = std::make_unique<bit_vector_stat>(
                    valid_edges_->convert_to<bit_vector_stat>()
                );
                break;
            }
            case Config::DYN: {
                valid_edges_ = std::make_unique<bit_vector_dyn>(
                    valid_edges_->convert_to<bit_vector_dyn>()
                );
                break;
            }
            case Config::SMALL: {
                valid_edges_ = std::make_unique<bit_vector_small>(
                    valid_edges_->convert_to<bit_vector_small>()
                );
                break;
            }
        }
    }

    boss_graph_->switch_state(new_state);
}

Config::StateType DBGSuccinct::get_state() const {
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::STAT
                || dynamic_cast<const bit_vector_stat*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::DYN
                || dynamic_cast<const bit_vector_dyn*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::SMALL
                || dynamic_cast<const bit_vector_small*>(valid_edges_.get()));

    return boss_graph_->get_state();
}

void DBGSuccinct::mask_dummy_kmers(size_t num_threads, bool with_pruning) {
    valid_edges_.reset();

    //TODO: use sdsl::bit_vector as mask
    std::vector<bool> vector_mask = with_pruning
        ? boss_graph_->prune_and_mark_all_dummy_edges(num_threads)
        : boss_graph_->mark_all_dummy_edges(num_threads);

    auto vector = to_sdsl(std::move(vector_mask));
    vector.flip();

    switch (get_state()) {
        case Config::STAT: {
            valid_edges_ = std::make_unique<bit_vector_stat>(std::move(vector));
            break;
        }
        case Config::DYN: {
            valid_edges_ = std::make_unique<bit_vector_dyn>(std::move(vector));
            break;
        }
        case Config::SMALL: {
            valid_edges_ = std::make_unique<bit_vector_small>(std::move(vector));
            break;
        }
    }

    assert(valid_edges_.get());
    assert(valid_edges_->size() == boss_graph_->num_edges() + 1);
    assert(!(*valid_edges_)[0]);
}

uint64_t DBGSuccinct::kmer_to_boss_index(node_index kmer_index) const {
    assert(kmer_index <= num_nodes());

    if (!valid_edges_.get() || !kmer_index)
        return kmer_index;

    return valid_edges_->select1(kmer_index);
}

DBGSuccinct::node_index DBGSuccinct::boss_to_kmer_index(uint64_t boss_index) const {
    assert(boss_index <= boss_graph_->num_edges());
    assert(!valid_edges_.get() || boss_index < valid_edges_->size());

    if (!valid_edges_.get() || !boss_index)
        return boss_index;

    if (!(*valid_edges_)[boss_index])
        return npos;

    return valid_edges_->rank1(boss_index);
}

bool DBGSuccinct::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const DBGSuccinct*>(&other)) {
        const auto& other_succ = *dynamic_cast<const DBGSuccinct*>(&other);
        if (boss_graph_.get() == other_succ.boss_graph_.get())
            return true;

        if (!boss_graph_.get() || !other_succ.boss_graph_.get())
            return false;

        return boss_graph_->equals_internally(*other_succ.boss_graph_, false);
    }

    throw std::runtime_error("Not implemented");
    return false;
}
