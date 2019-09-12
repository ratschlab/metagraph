#include "dbg_succinct.hpp"

#include <cassert>
#include <vector>
#include <algorithm>
#include <string>

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
    assert(in_graph(node));

    // dbg node is a boss edge
    BOSS::edge_index edge = boss_graph_->fwd(kmer_to_boss_index(node));
    return boss_to_kmer_index(
        boss_graph_->pick_edge(edge, boss_graph_->encode(next_char))
    );
}

// Traverse the incoming edge
node_index DBGSuccinct::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    // dbg node is a boss edge
    BOSS::edge_index edge = boss_graph_->bwd(kmer_to_boss_index(node));
    return boss_to_kmer_index(
        boss_graph_->pick_incoming_edge(edge, boss_graph_->encode(prev_char))
    );
}

template <class Callback>
inline void call_outgoing(const BOSS &boss,
                          uint64_t boss_edge,
                          const Callback &callback) {
    // no outgoing edges from the sink dummy nodes
    if (boss_edge > 1 && !boss.get_W(boss_edge))
        return;

    auto last = boss.fwd(boss_edge);
    auto first = boss.pred_last(last - 1) + 1;

    for (auto i = std::max(uint64_t(2), first); i <= last; ++i) {
        assert(boss.get_W(boss_edge) % boss.alph_size
                == boss.get_node_last_value(i));

        callback(i);
    }
}

void DBGSuccinct::call_outgoing_kmers(node_index node,
                                      const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    call_outgoing(*boss_graph_, kmer_to_boss_index(node), [&](auto i) {
        auto next = boss_to_kmer_index(i);
        if (next != npos)
            callback(next, boss_graph_->decode(boss_graph_->get_W(i)
                                % boss_graph_->alph_size));
    });
}

void DBGSuccinct::call_incoming_kmers(node_index node,
                                      const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));

    auto edge = kmer_to_boss_index(node);

    boss_graph_->call_incoming_to_target(boss_graph_->bwd(edge),
        [&](BOSS::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = boss_to_kmer_index(incoming_boss_edge);
            if (prev != npos)
                callback(prev, boss_graph_->get_minus_k_value(incoming_boss_edge,
                                                              get_k() - 2).first);
        }
    );
}

void DBGSuccinct::adjacent_outgoing_nodes(node_index node,
                                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing(*boss_graph_, kmer_to_boss_index(node), [&](auto i) {
        auto next = boss_to_kmer_index(i);
        if (next != npos)
            callback(next);
    });
}

void DBGSuccinct::adjacent_incoming_nodes(node_index node,
                                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    auto edge = kmer_to_boss_index(node);

    boss_graph_->call_incoming_to_target(boss_graph_->bwd(edge),
        [&](BOSS::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = boss_to_kmer_index(incoming_boss_edge);
            if (prev != npos)
                callback(prev);
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
            if (nodes_inserted)
                nodes_inserted->insert_bit(boss_to_kmer_index(i), true);
        }

    } else {
        boss_graph_->add_sequence(sequence, true);
    }

    assert(!valid_edges_.get() || !(*valid_edges_)[0]);
}

std::string DBGSuccinct::get_node_sequence(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    return boss_graph_->get_node_str(boss_edge)
            + boss_graph_->decode(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size);
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

template <class StringIt>
void DBGSuccinct
::call_nodes_with_suffix(StringIt begin,
                         StringIt end,
                         const std::function<void(node_index, uint64_t /* match length */)>& callback,
                         size_t min_match_length,
                         size_t max_num_allowed_matches) const {
    if (begin >= end || !max_num_allowed_matches)
        return;

    assert(get_k() >= static_cast<size_t>(end - begin));

    if (static_cast<size_t>(end - begin) < min_match_length)
        return;

    auto encoded = boss_graph_->encode(std::string(begin, end));
    auto index_range = boss_graph_->index_range(
        encoded.begin(),
        std::min(encoded.begin() + get_k() - 1, encoded.end())
    );

    if (std::get<0>(index_range) == 0 || std::get<1>(index_range) == 0)
        return;

    if (std::get<2>(index_range) == encoded.begin()) {
        callback(std::get<0>(index_range), 0);
        return;
    }

    // since we can only match up to get_k() - 1 in BOSS, check for this
    // case and simply pick the appropriate BOSS edge
    if (encoded.size() == get_k() && std::get<2>(index_range) + 1 == encoded.end()) {
        assert(std::get<0>(index_range) == std::get<1>(index_range));
        auto edge = boss_graph_->pick_edge(std::get<1>(index_range), encoded.back());
        if (edge) {
            auto kmer_index = boss_to_kmer_index(edge);
            if (kmer_index != npos) {
                callback(kmer_index, end - begin);
                return;
            }
        }
    }

    uint64_t match_size = std::get<2>(index_range) - encoded.begin();
    if (match_size < min_match_length)
        return;

    auto rank_first = boss_graph_->rank_last(std::get<0>(index_range));
    auto rank_last = boss_graph_->rank_last(std::get<1>(index_range));
    if (max_num_allowed_matches < std::numeric_limits<size_t>::max()) {
        std::vector<node_index> nodes;

        for (auto i = rank_first; i <= rank_last && nodes.size() <= max_num_allowed_matches; ++i) {
            boss_graph_->call_incoming_to_target(
                boss_graph_->bwd(boss_graph_->select_last(i)),
                [&](BOSS::edge_index incoming_edge_idx) {
                    auto kmer_index = boss_to_kmer_index(incoming_edge_idx);
                    if (kmer_index != npos)
                        nodes.emplace_back(kmer_index);
                }
            );
        }

        if (nodes.size() > max_num_allowed_matches)
            return;

        for (auto node : nodes) {
            callback(node, match_size);
        }
    } else {
        for (auto i = rank_first; i <= rank_last; ++i) {
            boss_graph_->call_incoming_to_target(
                boss_graph_->bwd(boss_graph_->select_last(i)),
                [&](BOSS::edge_index incoming_edge_idx) {
                    auto kmer_index = boss_to_kmer_index(incoming_edge_idx);
                    if (kmer_index != npos)
                        callback(kmer_index, match_size);
                }
            );
        }
    }
}

template void DBGSuccinct
::call_nodes_with_suffix<const char*>(const char*,
                                      const char*,
                                      const std::function<void(node_index, uint64_t)>&,
                                      size_t,
                                      size_t) const;

template void DBGSuccinct
::call_nodes_with_suffix<char*>(char*,
                                char*,
                                const std::function<void(node_index, uint64_t)>&,
                                size_t,
                                size_t) const;

template void DBGSuccinct
::call_nodes_with_suffix<std::string::const_iterator>(std::string::const_iterator,
                                                      std::string::const_iterator,
                                                      const std::function<void(node_index, uint64_t)>&,
                                                      size_t,
                                                      size_t) const;

template void DBGSuccinct
::call_nodes_with_suffix<std::string::iterator>(std::string::iterator,
                                                std::string::iterator,
                                                const std::function<void(node_index, uint64_t)>&,
                                                size_t,
                                                size_t) const;

void DBGSuccinct::traverse(node_index start,
                           const char *begin,
                           const char *end,
                           const std::function<void(node_index)> &callback,
                           const std::function<bool()> &terminate) const {
    assert(in_graph(start));
    assert(end >= begin);

    if (terminate())
        return;

    auto edge = kmer_to_boss_index(start);
    assert(edge);

    for (; begin != end && !terminate() && boss_graph_->get_W(edge); ++begin) {
        edge = boss_graph_->fwd(edge);
        edge = boss_graph_->pick_edge(edge, boss_graph_->encode(*begin));

        if (!edge)
            return;

        start = boss_to_kmer_index(edge);

        if (start != npos) {
            callback(start);
        } else {
            return;
        }
    }
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

void DBGSuccinct
::call_sequences(const std::function<void(const std::string&)> &callback) const {
    assert(boss_graph_.get());
    boss_graph_->call_sequences(callback);
}

void DBGSuccinct
::call_unitigs(const std::function<void(const std::string&)> &callback,
               size_t min_tip_size) const {
    assert(boss_graph_.get());
    boss_graph_->call_unitigs(callback, min_tip_size);
}

void DBGSuccinct
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    assert(boss_graph_.get());
    boss_graph_->call_kmers([&](auto index, const std::string &seq) {
        auto node = boss_to_kmer_index(index);
        assert(node != npos);
        callback(node, seq);
    });
}

void DBGSuccinct
::call_source_nodes(const std::function<void(node_index)> &callback) const {
    boss_graph_->call_start_edges([&](auto boss_edge) {
        auto node = boss_to_kmer_index(boss_edge);
        assert(node != npos);
        assert(!indegree(node));
        callback(node);
    });
}

size_t DBGSuccinct::outdegree(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) - 1;

    if (!boss_graph_->get_W(boss_edge)) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return 0;
    }

    auto last_target_kmer = boss_graph_->fwd(boss_edge);

    if (!boss_graph_->get_W(last_target_kmer)) {
        // There is a sink dummy target, hence this is the only outgoing edge
        // skip boss dummy sink edges
        return 0;
    }

    return last_target_kmer - boss_graph_->pred_last(last_target_kmer - 1);
}

bool DBGSuccinct::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) == 2;

    if (!boss_graph_->get_W(boss_edge)) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return false;
    }

    auto last_target_kmer = boss_graph_->fwd(boss_edge);

    if (!boss_graph_->get_W(last_target_kmer)) {
        // There is a sink dummy target, hence this is the only outgoing edge
        // skip boss dummy sink edges
        return false;
    }

    return boss_graph_->is_single_outgoing(last_target_kmer);
}

bool DBGSuccinct::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) > 2;

    if (!boss_graph_->get_W(boss_edge)) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return false;
    }

    return !boss_graph_->get_last(boss_graph_->fwd(boss_edge) - 1);
}

size_t DBGSuccinct::indegree(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return 1;

    auto x = boss_graph_->bwd(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    return boss_graph_->num_incoming_to_target(x) - !first_valid;
}

bool DBGSuccinct::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return false;

    auto x = boss_graph_->bwd(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    return !first_valid && boss_graph_->is_single_incoming(x);
}

bool DBGSuccinct::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = kmer_to_boss_index(node);

    if (boss_edge == 1)
        return false;

    auto x = boss_graph_->bwd(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    if (x + 1 == boss_graph_->get_W().size())
        return first_valid;

    if (first_valid)
        return boss_graph_->is_single_incoming(x);

    return boss_graph_->num_incoming_to_target(x) == 2;
}

uint64_t DBGSuccinct::num_nodes() const {
    return valid_edges_.get()
                ? valid_edges_->num_set_bits()
                : boss_graph_->num_edges();
}

bool DBGSuccinct::load_without_mask(const std::string &filename) {
    // release the old mask
    valid_edges_.reset();

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

    return true;
}

bool DBGSuccinct::load(const std::string &filename) {
    if (!load_without_mask(filename))
        return false;

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
        case Config::FAST: {
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
        || (boss_graph_->get_state() == Config::StateType::FAST
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
            case Config::FAST: {
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
                || boss_graph_->get_state() != Config::StateType::FAST
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

    auto vector_mask = with_pruning
        ? boss_graph_->prune_and_mark_all_dummy_edges(num_threads)
        : boss_graph_->mark_all_dummy_edges(num_threads);

    vector_mask.flip();

    switch (get_state()) {
        case Config::STAT: {
            valid_edges_ = std::make_unique<bit_vector_stat>(std::move(vector_mask));
            break;
        }
        case Config::FAST: {
            valid_edges_ = std::make_unique<bit_vector_stat>(std::move(vector_mask));
            break;
        }
        case Config::DYN: {
            valid_edges_ = std::make_unique<bit_vector_dyn>(std::move(vector_mask));
            break;
        }
        case Config::SMALL: {
            valid_edges_ = std::make_unique<bit_vector_small>(std::move(vector_mask));
            break;
        }
    }

    assert(valid_edges_.get());
    assert(valid_edges_->size() == boss_graph_->num_edges() + 1);
    assert(!(*valid_edges_)[0]);
}

void DBGSuccinct::reset_mask() {
    valid_edges_.reset();
}

uint64_t DBGSuccinct::kmer_to_boss_index(node_index kmer_index) const {
    assert(in_graph(kmer_index));

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
    if (get_k() != other.get_k()
            || num_nodes() != other.num_nodes()
            || is_canonical_mode() != other.is_canonical_mode())
        return false;

    if (dynamic_cast<const DBGSuccinct*>(&other)) {
        const auto &other_succ = *dynamic_cast<const DBGSuccinct*>(&other);

        if (this == &other_succ)
            return true;

        // only one of the mask vectors is defined
        if (bool(valid_edges_.get()) != bool(other_succ.valid_edges_.get()))
            return false;

        // different mask vectors
        if (valid_edges_.get() && *valid_edges_ != *other_succ.valid_edges_)
            return false;

        // TODO: what if graphs have same real nodes, different
        // sets of dummy nodes but they all are masked out?

        assert(boss_graph_.get() && other_succ.boss_graph_.get());

        return boss_graph_->equals_internally(*other_succ.boss_graph_, false);
    }

    throw std::runtime_error("Not implemented");
    return false;
}


const std::string& DBGSuccinct::alphabet() const {
    return boss_graph_->alphabet;
}

void DBGSuccinct::print(std::ostream &out) const {
    auto vertex_header = std::string("Vertex");
    vertex_header.resize(get_k() - 1, ' ');

    out << "BOSS" << "\t" << "L"
                  << "\t" << vertex_header
                  << "\t" << "W";

    if (valid_edges_.get())
        out << "\t" << "Index" << "\t" << "Valid";

    out << std::endl;

    const auto &boss = get_boss();

    uint64_t valid_count = 0;

    for (uint64_t i = 1; i <= boss.num_edges(); i++) {
        out << i << "\t" << boss.get_last(i)
                 << "\t" << boss.get_node_str(i)
                 << "\t" << boss.decode(boss.get_W(i) % boss.alph_size)
                         << (boss.get_W(i) >= boss.alph_size
                                 ? "-"
                                 : "");

        if (valid_edges_.get()) {
            bool valid = (*valid_edges_)[i];
            valid_count += valid;
            out << "\t" << (valid ? valid_count : 0)
                << "\t" << valid;
        }

        out << std::endl;
    }
}

bool DBGSuccinct::in_graph(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    std::ignore = node;
    return true;
}
