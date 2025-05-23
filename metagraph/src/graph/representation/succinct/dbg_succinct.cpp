#include "dbg_succinct.hpp"

#include <cassert>
#include <vector>
#include <algorithm>
#include <string>
#include <filesystem>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/serialization.hpp"
#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"


namespace mtg {
namespace graph {

using namespace mtg::graph::boss;

using mtg::common::logger;

typedef DBGSuccinct::node_index node_index;


DBGSuccinct::DBGSuccinct(size_t k, Mode mode)
      : boss_graph_(std::make_unique<BOSS>(k - 1)),
        mode_(mode) {
    if (k < 2) {
        logger->error("For succinct graph, k must be at least 2");
        exit(1);
    }
}

DBGSuccinct::DBGSuccinct(BOSS *boss_graph, Mode mode)
      : boss_graph_(boss_graph),
        mode_(mode) {}

size_t DBGSuccinct::get_k() const {
    return boss_graph_->get_k() + 1;
}

// Check whether graph contains fraction of nodes from the sequence
bool DBGSuccinct::find(std::string_view sequence,
                       double discovery_fraction) const {
    if (sequence.length() < get_k())
        return false;

    const size_t num_kmers = sequence.length() - get_k() + 1;
    const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    auto is_invalid = get_missing_kmer_skipper(bloom_filter_.get(), sequence);

    boss_graph_->map_to_edges(sequence,
        [&](auto edge) {
            if (edge) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; },
        [&]() {
            if (!is_invalid())
                return false;

            num_kmers_missing++;
            return true;
        }
    );

    return num_kmers_missing <= max_kmers_missing;
}

// Traverse the outgoing edge
node_index DBGSuccinct::traverse(node_index node, char next_char) const {
    assert(in_graph(node));

    // return npos if the character is invalid
    if (boss_graph_->encode(next_char) == boss_graph_->alph_size)
        return npos;

    // dbg node is a boss edge
    BOSS::edge_index boss_edge = node;
    boss_edge = boss_graph_->fwd(boss_edge);
    return validate_edge(
        boss_graph_->pick_edge(boss_edge, boss_graph_->encode(next_char))
    );
}

// Traverse the incoming edge
node_index DBGSuccinct::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    // dbg node is a boss edge
    BOSS::edge_index edge = boss_graph_->bwd(node);
    return validate_edge(
        boss_graph_->pick_incoming_edge(edge, boss_graph_->encode(prev_char))
    );
}

template <class Callback>
inline void call_outgoing(const BOSS &boss,
                          uint64_t boss_edge,
                          const Callback &callback) {
    // no outgoing edges from the sink dummy nodes
    BOSS::TAlphabet w = 0;
    if (boss_edge > 1 && !(w = boss.get_W(boss_edge)))
        return;

    auto last = boss.fwd(boss_edge, w % boss.alph_size);
    auto first = boss.pred_last(last - 1) + 1;

    for (auto i = std::max(uint64_t(2), first); i <= last; ++i) {
        assert(w % boss.alph_size == boss.get_node_last_value(i));

        callback(i);
    }
}

void DBGSuccinct::call_outgoing_kmers(node_index node,
                                      const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    call_outgoing(*boss_graph_, node, [&](auto i) {
        auto next = i;
        if (in_graph(next))
            callback(next, boss_graph_->decode(boss_graph_->get_W(i)
                                % boss_graph_->alph_size));
    });
}

void DBGSuccinct::call_incoming_kmers(node_index node,
                                      const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));

    auto edge = node;

    boss_graph_->call_incoming_to_target(boss_graph_->bwd(edge),
        boss_graph_->get_node_last_value(edge),
        [&](BOSS::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = incoming_boss_edge;
            if (in_graph(prev)) {
                callback(prev,
                    boss_graph_->decode(
                        boss_graph_->get_minus_k_value(incoming_boss_edge, get_k() - 2).first
                    )
                );
            }
        }
    );
}

void DBGSuccinct::adjacent_outgoing_nodes(node_index node,
                                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing(*boss_graph_, node, [&](auto i) {
        auto next = i;
        if (in_graph(next))
            callback(next);
    });
}

void DBGSuccinct::adjacent_incoming_nodes(node_index node,
                                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    auto edge = node;

    boss_graph_->call_incoming_to_target(boss_graph_->bwd(edge),
        boss_graph_->get_node_last_value(edge),
        [&](BOSS::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = incoming_boss_edge;
            if (in_graph(prev))
                callback(prev);
        }
    );
}

void DBGSuccinct::call_nodes(const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate,
                             size_t num_threads,
                             size_t batch_size) const {
    if (valid_edges_) {
        size_t block_size = max_index() / num_threads;

        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (size_t begin = 1; begin <= max_index(); begin += block_size) {
            if (terminate())
                continue;

            size_t end = std::min(begin + block_size,
                                  static_cast<size_t>(max_index() + 1));
            try {
                valid_edges_->call_ones_in_range(begin, end, [&](uint64_t i) {
                    callback(i);
                    if (terminate())
                        throw early_term();
                });
            } catch (early_term&) {}
        }
    } else if (!terminate()) {
        try {
            call_sequences([&](const std::string&, const auto &path) {
                for (node_index node : path) {
                    if (terminate())
                        throw early_term();

                    callback(node);
                }
            }, num_threads);
        } catch (early_term&) {}
    }
}

void DBGSuccinct::add_sequence(std::string_view sequence,
                               const std::function<void(node_index)> &on_insertion) {
    if (sequence.size() < get_k())
        return;

    std::vector<uint64_t> boss_edges_inserted;
    boss_edges_inserted.reserve((sequence.size() - get_k() + 1) * 2);

    // insert forward sequence
    boss_graph_->add_sequence(sequence, true, &boss_edges_inserted);

    if (mode_ == CANONICAL) {
        // insert reverse complement sequence as well,
        // to have all canonical k-mers in graph
        std::string rev_compl(sequence.begin(), sequence.end());
        reverse_complement(rev_compl.begin(), rev_compl.end());

        boss_graph_->add_sequence(rev_compl, true, &boss_edges_inserted);
    }

    for (uint64_t new_boss_edge : boss_edges_inserted) {
        if (get_state() != BOSS::State::DYN)
            throw std::runtime_error("representation must be dynamic");

        // update bitmask with valid k-mers -- assume all inserted k-mers valid
        // TODO: detect dummy BOSS edges and insert zeros to mask them out
        if (valid_edges_)
            dynamic_cast<bit_vector_dyn&>(*valid_edges_).insert_bit(new_boss_edge, 1);

        // Call all new nodes inserted including the dummy ones, unless they
        // are masked out.
        on_insertion(validate_edge(new_boss_edge));
    }

    assert(!valid_edges_.get() || !(*valid_edges_)[0]);
    assert(!valid_edges_.get() || valid_edges_->size() == boss_graph_->num_edges() + 1);

    if (bloom_filter_)
        bloom_filter_->add_sequence(sequence);
}

std::string DBGSuccinct::get_node_sequence(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    return boss_graph_->get_node_str(boss_edge)
            + boss_graph_->decode(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size);
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are not mapped to canonical ones
void DBGSuccinct::map_to_nodes_sequentially(std::string_view sequence,
                                            const std::function<void(node_index)> &callback,
                                            const std::function<bool()> &terminate) const {
    if (sequence.size() < get_k())
        return;

    auto is_missing = get_missing_kmer_skipper(bloom_filter_.get(), sequence);

    boss_graph_->map_to_edges(
        sequence,
        [&](BOSS::edge_index i) { callback(validate_edge(i)); },
        terminate,
        [&]() {
            if (!is_missing())
                return false;

            callback(npos);
            return true;
        }
    );
}

void DBGSuccinct
::call_nodes_with_suffix_matching_longest_prefix(
            std::string_view str,
            std::function<void(node_index, uint64_t /* match length */)> callback,
            size_t min_match_length,
            size_t max_num_allowed_matches) const {
    assert(str.size() <= get_k());
    if (!max_num_allowed_matches || str.size() < min_match_length)
        return;

    auto encoded = boss_graph_->encode(str);

    if (std::find(encoded.begin(), encoded.end(),
                  boss_graph_->alph_size) != encoded.end()) {
        return;
    }

    auto [first, last, end] = boss_graph_->index_range(
        encoded.begin(),
        std::min(encoded.begin() + get_k() - 1, encoded.end())
    );
    size_t match_size = end - encoded.begin();

    // since we can only match up to get_k() - 1 in BOSS, check for this
    // case and simply pick the appropriate BOSS edge
    if (str.size() == get_k() && match_size + 1 == get_k()) {
        assert(first == last);
        auto edge = boss_graph_->pick_edge(last, encoded.back());
        if (edge) {
            auto kmer_index = edge;
            if (in_graph(kmer_index)) {
                assert(str.size() == get_k());
                assert(get_node_sequence(kmer_index) == str);
                callback(kmer_index, get_k());
                return;
            }
        }
    }

    if (match_size < min_match_length)
        return;

    auto rank_first = boss_graph_->rank_last(first);
    auto rank_last = boss_graph_->rank_last(last);
    // TODO: rewrite this, call first N nodes and discard the large batches
    // in the caller (if needed at all)
    if (max_num_allowed_matches < std::numeric_limits<size_t>::max()) {
        std::vector<node_index> nodes;

        for (auto i = rank_first; i <= rank_last && nodes.size() <= max_num_allowed_matches; ++i) {
            BOSS::edge_index e = boss_graph_->select_last(i);
            boss_graph_->call_incoming_to_target(boss_graph_->bwd(e),
                boss_graph_->get_node_last_value(e),
                [&](BOSS::edge_index incoming_edge_idx) {
                    auto kmer_index = incoming_edge_idx;
                    if (in_graph(kmer_index)) {
                        assert(get_node_sequence(kmer_index).substr(get_k() - match_size)
                            == str.substr(0, match_size));
                        nodes.emplace_back(kmer_index);
                    }
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
            BOSS::edge_index e = boss_graph_->select_last(i);
            boss_graph_->call_incoming_to_target(boss_graph_->bwd(e),
                boss_graph_->get_node_last_value(e),
                [&](BOSS::edge_index incoming_edge_idx) {
                    auto kmer_index = incoming_edge_idx;
                    if (in_graph(kmer_index)) {
                        assert(get_node_sequence(kmer_index).substr(get_k() - match_size)
                            == str.substr(0, match_size));
                        callback(kmer_index, match_size);
                    }
                }
            );
        }
    }
}

void DBGSuccinct::traverse(node_index start,
                           const char *begin,
                           const char *end,
                           const std::function<void(node_index)> &callback,
                           const std::function<bool()> &terminate) const {
    assert(in_graph(start));
    assert(end >= begin);

    if (terminate())
        return;

    auto edge = start;
    assert(edge);

    BOSS::TAlphabet w;
    for (; begin != end && !terminate() && (w = boss_graph_->get_W(edge)); ++begin) {
        // stop traversal if the character is invalid
        if (boss_graph_->encode(*begin) == boss_graph_->alph_size)
            return;

        edge = boss_graph_->fwd(edge, w % boss_graph_->alph_size);
        edge = boss_graph_->pick_edge(edge, boss_graph_->encode(*begin));

        start = edge;
        if (!in_graph(start))
            return;

        callback(start);
    }
}

// Map sequence k-mers to the canonical graph nodes
// and run callback for each node until the termination condition is satisfied
void DBGSuccinct::map_to_nodes(std::string_view sequence,
                               const std::function<void(node_index)> &callback,
                               const std::function<bool()> &terminate) const {
    if (sequence.size() < get_k())
        return;

    auto is_missing = get_missing_kmer_skipper(bloom_filter_.get(), sequence);

    if (mode_ == CANONICAL) {
        std::string sequence_rev_compl(sequence.begin(), sequence.end());
        reverse_complement(sequence_rev_compl.begin(), sequence_rev_compl.end());

        std::vector<BOSS::edge_index> boss_edges;
        boss_edges.resize(sequence.size() - get_k() + 1);

        auto jt = boss_edges.begin();
        boss_graph_->map_to_edges(sequence,
            [&](auto index) { *jt = index; ++jt; },
            []() { return false; },
            [&]() {
                if (!is_missing())
                    return false;

                ++jt;
                return true;
            }
        );

        assert(jt == boss_edges.end());

        auto it = boss_edges.rbegin();
        boss_graph_->map_to_edges(sequence_rev_compl,
            [&](auto rc_index) {
                assert(it < boss_edges.rend());
                *it = std::min(*it, rc_index);
                ++it;
            },
            []() { return false; },
            [&]() {
                assert(it < boss_edges.rend());
                // if a k-mer is missing, skip its reverse compliment, as it's missing too.
                if (!*it) {
                    ++it;
                    return true;
                } else {
                    return false;
                }
            }
        );

        for (size_t i = 0; i < boss_edges.size() && !terminate(); ++i) {
            // the definition of a canonical k-mer is redefined:
            //      use k-mer with smaller index in the BOSS table.
            callback(validate_edge(boss_edges[i]));
        }

    } else {
        boss_graph_->map_to_edges(
            sequence,
            [&](BOSS::edge_index i) { callback(validate_edge(i)); },
            terminate,
            [&]() {
                if (!is_missing())
                    return false;

                callback(npos);
                return true;
            }
        );
    }
}

void DBGSuccinct::call_sequences(const CallPath &callback,
                                 size_t num_threads,
                                 bool kmers_in_single_form,
                                 bool verbose) const {
    assert(boss_graph_.get());
    boss_graph_->call_sequences(
        [&](std::string&& seq, auto&& path) {
            for (auto &node : path) {
                node = validate_edge(node);
            }
            callback(std::move(seq), std::move(path));
        },
        num_threads,
        kmers_in_single_form,
        verbose
    );
}

void DBGSuccinct::call_unitigs(const CallPath &callback,
                               size_t num_threads,
                               size_t min_tip_size,
                               bool kmers_in_single_form) const {
    assert(boss_graph_.get());
    boss_graph_->call_unitigs(
        [&](std::string&& seq, auto&& path) {
            for (auto &node : path) {
                node = validate_edge(node);
            }
            callback(std::move(seq), std::move(path));
        },
        num_threads,
        min_tip_size,
        kmers_in_single_form
    );
}

void DBGSuccinct
::call_kmers(const std::function<void(node_index, const std::string&)> &callback,
             const std::function<bool()> &stop_early) const {
    assert(boss_graph_.get());
    boss_graph_->call_kmers([&](auto index, const std::string &seq) {
        auto node = index;
        assert(in_graph(node));
        callback(node, seq);
    }, stop_early);
}

void DBGSuccinct
::call_source_nodes(const std::function<void(node_index)> &callback) const {
    boss_graph_->call_start_edges([&](auto boss_edge) {
        auto node = boss_edge;
        assert(in_graph(node));
        assert(!indegree(node));
        callback(node);
    });
}

size_t DBGSuccinct::outdegree(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) - 1;

    BOSS::TAlphabet d = boss_graph_->get_W(boss_edge) % boss_graph_->alph_size;
    if (!d) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return 0;
    }

    auto last_target_kmer = boss_graph_->fwd(boss_edge, d);

    if (!(valid_edges_ ? (*valid_edges_)[last_target_kmer]
                       : boss_graph_->get_W(last_target_kmer))) {
        // There is a sink dummy target, hence this is the only outgoing edge
        // skip boss dummy sink edges
        return 0;
    }

    return last_target_kmer - boss_graph_->pred_last(last_target_kmer - 1);
}

bool DBGSuccinct::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) == 2;

    BOSS::TAlphabet d = boss_graph_->get_W(boss_edge) % boss_graph_->alph_size;
    if (!d) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return false;
    }

    auto last_target_kmer = boss_graph_->fwd(boss_edge, d);

    if (!(valid_edges_ ? (*valid_edges_)[last_target_kmer]
                       : boss_graph_->get_W(last_target_kmer))) {
        // There is a sink dummy target, hence this is the only outgoing edge
        // skip boss dummy sink edges
        return false;
    }

    return boss_graph_->is_single_outgoing(last_target_kmer);
}

bool DBGSuccinct::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return boss_graph_->succ_last(1) > 2;

    BOSS::TAlphabet d = boss_graph_->get_W(boss_edge) % boss_graph_->alph_size;
    if (!d) {
        // |node| is a sink dummy boss edge, hence has no outgoing edges
        return false;
    }

    return !boss_graph_->get_last(boss_graph_->fwd(boss_edge, d) - 1);
}

size_t DBGSuccinct::indegree(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return 1;

    auto x = boss_graph_->bwd(boss_edge);
    BOSS::TAlphabet w = boss_graph_->get_node_last_value(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    return boss_graph_->num_incoming_to_target(x, w) - !first_valid;
}

bool DBGSuccinct::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return false;

    auto x = boss_graph_->bwd(boss_edge);
    BOSS::TAlphabet w = boss_graph_->get_node_last_value(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    return !first_valid && boss_graph_->is_single_incoming(x, w);
}

bool DBGSuccinct::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    auto boss_edge = node;

    if (boss_edge == 1)
        return false;

    auto x = boss_graph_->bwd(boss_edge);
    BOSS::TAlphabet w = boss_graph_->get_node_last_value(boss_edge);

    size_t first_valid = !valid_edges_.get() || (*valid_edges_)[x];

    if (x + 1 == boss_graph_->get_W().size())
        return first_valid;

    if (first_valid)
        return boss_graph_->is_single_incoming(x, w);

    return boss_graph_->num_incoming_to_target(x, w) == 2;
}

uint64_t DBGSuccinct::num_nodes() const {
    return valid_edges_.get()
                ? valid_edges_->num_set_bits()
                : boss_graph_->num_edges();
}

uint64_t DBGSuccinct::max_index() const {
    return boss_graph_->num_edges();
}

bool DBGSuccinct::load_without_mask(const std::string &filename) {
    // release the old mask
    valid_edges_.reset();

    {
        std::unique_ptr<std::ifstream> in
            = utils::open_ifstream(utils::make_suffix(filename, kExtension));

        if (!boss_graph_->load(*in))
            return false;

        mode_ = static_cast<Mode>(load_number(*in));

        if (!boss_graph_->load_suffix_ranges(*in))
            logger->warn("No index for node ranges could be loaded");
    }

    return true;
}

bool DBGSuccinct::load(const std::string &filename) {
    if (!load_without_mask(filename))
        return false;

    auto prefix = utils::remove_suffix(filename, kExtension);

    std::unique_ptr<std::ifstream> in = utils::open_ifstream(prefix + kDummyMaskExtension);
    if (!in->good())
        return true;

    // initialize a new vector
    switch (get_state()) {
        case BOSS::State::STAT: {
            valid_edges_.reset(new bit_vector_small());
            break;
        }
        case BOSS::State::FAST: {
            valid_edges_.reset(new bit_vector_stat());
            break;
        }
        case BOSS::State::DYN: {
            valid_edges_.reset(new bit_vector_dyn());
            break;
        }
        case BOSS::State::SMALL: {
            valid_edges_.reset(new bit_vector_small());
            break;
        }
    }

    // load the mask of valid edges (all non-dummy including npos 0)
    if (!valid_edges_->load(*in)) {
        std::cerr << "Error: Can't load dummy edge mask." << std::endl;
        return false;
    }

    if (valid_edges_->size() != boss_graph_->num_edges() + 1 || (*valid_edges_)[0]) {
        std::cerr << "Error: Edge mask is not compatible with graph." << std::endl;
        return false;
    }

    if (std::filesystem::exists(prefix + kBloomFilterExtension)) {
        std::unique_ptr<std::ifstream> bloom_in = utils::open_ifstream(prefix + kBloomFilterExtension);
        if (!bloom_filter_)
            bloom_filter_ = std::make_unique<kmer::KmerBloomFilter<>>(get_k(), mode_ == CANONICAL);

        if (!bloom_filter_->load(*bloom_in)) {
            std::cerr << "Error: failed to load Bloom filter from " + prefix + kBloomFilterExtension << std::endl;
            return false;
        }

        assert(bloom_filter_);

        if (bloom_filter_->is_canonical_mode() != (get_mode() == CANONICAL)) {
            std::cerr << "Error: Bloom filter and graph in incompatible modes" << std::endl
                      << "Bloom filter: " << (bloom_filter_->is_canonical_mode() ? "not " : "") << "canonical" << std::endl
                      << "Graph: " << static_cast<int>(get_mode()) << std::endl;
            return false;
        }

        if (bloom_filter_->get_k() != get_k()) {
            std::cerr << "Error: mismatched k between Bloom filter and graph" << std::endl
                      << "Bloom filter: " << bloom_filter_->get_k() << std::endl
                      << "Graph: " << get_k() << std::endl;
            return false;
        }

        logger->trace("Bloom filter loaded from {} successfully",
                      prefix + kBloomFilterExtension);
    }

    return true;
}

void DBGSuccinct::serialize(const std::string &filename) const {
    const std::string &prefix = utils::remove_suffix(filename, kExtension);

    // Clear any existing Bloom filters
    std::filesystem::remove(prefix + kBloomFilterExtension);

    {
        const std::string out_filename = prefix + kExtension;
        std::ofstream out = utils::open_new_ofstream(out_filename);
        boss_graph_->serialize(out);
        serialize_number(out, static_cast<int>(mode_));

        boss_graph_->serialize_suffix_ranges(out);

        if (!out.good())
            throw std::ios_base::failure("Can't write to file " + out_filename);
    }

    auto serialize_valid_edges = [&](auto &&valid_edges) {
        assert((boss_graph_->get_state() == BOSS::State::STAT
                    && dynamic_cast<const bit_vector_small*>(valid_edges.get()))
            || (boss_graph_->get_state() == BOSS::State::FAST
                    && dynamic_cast<const bit_vector_stat*>(valid_edges.get()))
            || (boss_graph_->get_state() == BOSS::State::DYN
                    && dynamic_cast<const bit_vector_dyn*>(valid_edges.get()))
            || (boss_graph_->get_state() == BOSS::State::SMALL
                    && dynamic_cast<const bit_vector_small*>(valid_edges.get())));

        const auto out_filename = prefix + kDummyMaskExtension;
        std::ofstream out = utils::open_new_ofstream(out_filename);
        if (!out.good())
            throw std::ios_base::failure("Can't write to file " + out_filename);

        valid_edges->serialize(out);
    };

    if (valid_edges_)
        serialize_valid_edges(valid_edges_);

    if (bloom_filter_) {
        std::ofstream bloom_out = utils::open_new_ofstream(prefix + kBloomFilterExtension);
        if (!bloom_out.good())
            throw std::ios_base::failure("Can't write to file " + prefix + kBloomFilterExtension);

        bloom_filter_->serialize(bloom_out);
    }
}

void DBGSuccinct::serialize(boss::BOSS::Chunk&& chunk,
                            const std::string &filename,
                            Mode mode, BOSS::State state) {
    const std::string &prefix = utils::remove_suffix(filename, kExtension);

    std::filesystem::remove(prefix + kBloomFilterExtension);
    std::filesystem::remove(prefix + kDummyMaskExtension);

    const std::string &fname = prefix + kExtension;
    std::ofstream out = utils::open_new_ofstream(fname);
    boss::BOSS::serialize(std::move(chunk), out, state);
    serialize_number(out, static_cast<int>(mode));
    serialize_number(out, 0); // suffix ranges are not indexed
    if (!out.good())
        throw std::ios_base::failure("Can't write to file " + fname);
}

void DBGSuccinct::switch_state(BOSS::State new_state) {
    if (get_state() == new_state)
        return;

    if (valid_edges_.get()) {
        switch (new_state) {
            case BOSS::State::STAT: {
                valid_edges_ = std::make_unique<bit_vector_small>(
                    valid_edges_->convert_to<bit_vector_small>()
                );
                break;
            }
            case BOSS::State::FAST: {
                valid_edges_ = std::make_unique<bit_vector_stat>(
                    valid_edges_->convert_to<bit_vector_stat>()
                );
                break;
            }
            case BOSS::State::DYN: {
                valid_edges_ = std::make_unique<bit_vector_dyn>(
                    valid_edges_->convert_to<bit_vector_dyn>()
                );
                break;
            }
            case BOSS::State::SMALL: {
                valid_edges_ = std::make_unique<bit_vector_small>(
                    valid_edges_->convert_to<bit_vector_small>()
                );
                break;
            }
        }
    }

    boss_graph_->switch_state(new_state);
}

BOSS::State DBGSuccinct::get_state() const {
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != BOSS::State::STAT
                || dynamic_cast<const bit_vector_small*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != BOSS::State::FAST
                || dynamic_cast<const bit_vector_stat*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != BOSS::State::DYN
                || dynamic_cast<const bit_vector_dyn*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != BOSS::State::SMALL
                || dynamic_cast<const bit_vector_small*>(valid_edges_.get()));

    return boss_graph_->get_state();
}

std::unique_ptr<bit_vector> DBGSuccinct::generate_valid_kmer_mask(size_t num_threads, bool with_pruning) const {
    auto vector_mask = with_pruning
        ? boss_graph_->prune_and_mark_all_dummy_edges(num_threads)
        : boss_graph_->mark_all_dummy_edges(num_threads);

    vector_mask.flip();

    switch (get_state()) {
        case BOSS::State::STAT:
            return std::make_unique<bit_vector_small>(std::move(vector_mask));
        case BOSS::State::FAST:
            return std::make_unique<bit_vector_stat>(std::move(vector_mask));
        case BOSS::State::DYN:
            return std::make_unique<bit_vector_dyn>(std::move(vector_mask));
        case BOSS::State::SMALL:
            return std::make_unique<bit_vector_small>(std::move(vector_mask));
        default:
            throw std::runtime_error("Invalid state");
    }
}

void DBGSuccinct::mask_dummy_kmers(size_t num_threads, bool with_pruning) {
    valid_edges_.reset();

    valid_edges_ = generate_valid_kmer_mask(num_threads, with_pruning);

    assert(valid_edges_.get());
    assert(valid_edges_->size() == boss_graph_->num_edges() + 1);
    assert(!(*valid_edges_)[0]);
}

bool DBGSuccinct::in_graph(node_index node) const {
    return DeBruijnGraph::in_graph(node) && (!valid_edges_ || (*valid_edges_)[node]);
}
node_index DBGSuccinct::validate_edge(node_index node) const {
    return in_graph(node) ? node : npos;
}
node_index DBGSuccinct::select_node(uint64_t rank) const {
    assert(rank <= num_nodes());

    if (!valid_edges_.get() || !rank)
        return rank;

    return valid_edges_->select1(rank);
}

uint64_t DBGSuccinct::rank_node(node_index node) const {
    assert(node <= max_index());

    if (!valid_edges_.get() || !node)
        return node;

    if (!(*valid_edges_)[node])
        return npos;

    return valid_edges_->rank1(node);
}

void DBGSuccinct
::initialize_bloom_filter_from_fpr(double false_positive_rate,
                                   uint32_t max_num_hash_functions) {
    bloom_filter_ = std::make_unique<kmer::KmerBloomFilter<>>(
        get_k(),
        mode_ == CANONICAL,
        BloomFilter::optim_size(false_positive_rate, num_nodes()),
        num_nodes(),
        std::min(max_num_hash_functions, BloomFilter::optim_h(false_positive_rate))
    );

    std::mutex seq_mutex;
    bloom_filter_->add_sequences([&](const auto &callback) {
        call_sequences(
            [&](const auto &sequence, const auto &) {
                std::lock_guard<std::mutex> lock(seq_mutex);
                callback(sequence);
            },
            get_num_threads(),
            mode_ == CANONICAL
        );
    });
}

void DBGSuccinct
::initialize_bloom_filter(double bits_per_kmer,
                          uint32_t max_num_hash_functions) {
    bloom_filter_ = std::make_unique<kmer::KmerBloomFilter<>>(
        get_k(),
        mode_ == CANONICAL,
        bits_per_kmer * num_nodes(),
        num_nodes(),
        max_num_hash_functions
    );

    std::mutex seq_mutex;
    bloom_filter_->add_sequences([&](const auto &callback) {
        call_sequences(
            [&](const auto &sequence, const auto &) {
                std::lock_guard<std::mutex> lock(seq_mutex);
                callback(sequence);
            },
            get_num_threads(),
            mode_ == CANONICAL
        );
    });
}

bool DBGSuccinct::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || num_nodes() != other.num_nodes()
            || get_mode() != other.get_mode())
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
    std::string vertex_header = "Vertex";
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
        BOSS::TAlphabet w = boss.get_W(i);
        assert(w != boss.alph_size);
        out << i << "\t" << boss.get_last(i)
                 << "\t" << boss.get_node_str(i)
                 << "\t" << boss.decode(w % boss.alph_size)
                         << (w > boss.alph_size ? "-" : "");

        if (valid_edges_.get()) {
            bool valid = (*valid_edges_)[i];
            valid_count += valid;
            out << "\t" << (valid ? valid_count : 0)
                << "\t" << valid;
        }

        out << std::endl;
    }
}

} // namespace graph
} // namespace mtg
