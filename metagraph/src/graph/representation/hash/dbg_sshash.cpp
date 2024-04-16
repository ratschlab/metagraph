#include "dbg_sshash.hpp"

#include <query/streaming_query_regular_parsing.hpp>


namespace mtg {
namespace graph {

constexpr DeBruijnGraph::node_index sshash_to_graph_index(uint64_t idx) {
    return idx + 1;
}
constexpr uint64_t graph_index_to_sshash(DeBruijnGraph::node_index idx) {
    return idx - 1;
}

DBGSSHash::DBGSSHash(size_t k) : k_(k) {}

DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k, Mode mode)
    : k_(k), mode_(mode) {
    sshash::build_configuration build_config;
    build_config.k = k;
    // quick fix for value of m... k/2 but odd
    build_config.m = (k_ + 1) / 2;
    if (build_config.m % 2 == 0)
        build_config.m++;
    dict_.build(input_filename, build_config);
}

std::string DBGSSHash::file_extension() const {
    return kExtension;
}
size_t DBGSSHash::get_k() const {
    return k_;
}
DeBruijnGraph::Mode DBGSSHash::get_mode() const {
    return mode_;
}

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)>& on_insertion) {
    throw std::runtime_error("adding sequences not implemented");
}

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)>& callback,
                             const std::function<bool()>& terminate) const {
    if (mode_ != CANONICAL) {
        map_to_nodes_sequentially(sequence, callback, terminate);
    } else {
        map_to_nodes_with_rc(
                sequence, [&](node_index node, bool) { callback(node); }, terminate);
    }
}

void DBGSSHash ::map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)>& callback,
                                           const std::function<bool()>& terminate) const {
    if (terminate() || sequence.size() < k_)
        return;

    kmer_t uint_kmer = sshash::util::string_to_uint_kmer<kmer_t>(sequence.data(), k_ - 1);
    uint_kmer.pad_char();
    for (size_t i = k_ - 1; i < sequence.size() && !terminate(); ++i) {
        uint_kmer.drop_char();
        uint_kmer.kth_char_or(k_ - 1, kmer_t::char_to_uint(sequence[i]));
        callback(sshash_to_graph_index(dict_.lookup_uint(uint_kmer, false)));
    }
}

void DBGSSHash ::map_to_nodes_with_rc(std::string_view sequence,
                                      const std::function<void(node_index, bool)>& callback,
                                      const std::function<bool()>& terminate) const {
    sshash::streaming_query_regular_parsing<kmer_t> streamer(&dict_);
    streamer.start();
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
        const char* kmer = sequence.data() + i;
        auto res = streamer.lookup_advanced(kmer);
        callback(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    std::string string_kmer = DBGSSHash::get_node_sequence(node);
    kmer_t new_kmer = sshash::util::string_to_uint_kmer<kmer_t>(string_kmer.c_str(), k_);
    new_kmer.drop_char();
    new_kmer.kth_char_or(k_ - 1, kmer_t::char_to_uint(next_char));
    uint64_t sshash_id
            = dict_.lookup_advanced_uint(new_kmer, /*check_reverse_complement*/ false).kmer_id;
    return sshash_to_graph_index(sshash_id);
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    std::string string_kmer = DBGSSHash::get_node_sequence(node);
    kmer_t new_kmer = sshash::util::string_to_uint_kmer<kmer_t>(string_kmer.c_str(), k_);
    new_kmer.append_char(kmer_t::char_to_uint(prev_char));
    new_kmer.take_chars(k_);
    uint64_t sshash_id
            = dict_.lookup_advanced_uint(new_kmer, /*check_reverse_complement*/ false).kmer_id;
    return sshash_to_graph_index(sshash_id);
}

void DBGSSHash ::adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

void DBGSSHash ::adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGSSHash ::call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_forward_neighbours(kmer.c_str(), false);
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
            callback(sshash_to_graph_index(nb.forward[i].kmer_id), kmer_t::alphabet[i]);
        }
    }
}

void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_backward_neighbours(kmer.c_str(), false);
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
            callback(sshash_to_graph_index(nb.backward[i].kmer_id), kmer_t::alphabet[i]);
        }
    }
}

void DBGSSHash ::call_outgoing_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_forward_neighbours(kmer.c_str(), false);
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
            callback(sshash_to_graph_index(nb.forward[i].kmer_id), kmer_t::alphabet[i],
                     nb.forward[i].kmer_orientation);
        }
    }
}


void DBGSSHash ::call_incoming_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_backward_neighbours(kmer.c_str(), false);
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
            callback(sshash_to_graph_index(nb.backward[i].kmer_id), kmer_t::alphabet[i],
                     nb.backward[i].kmer_orientation);
        }
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_forward_neighbours(kmer.c_str(), false);
    return std::count_if(begin(nb.forward), end(nb.forward), [](auto const& lookup_result) {
        return lookup_result.kmer_id != sshash::constants::invalid_uint64;
    });
}

size_t DBGSSHash::indegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood<kmer_t> nb = dict_.kmer_backward_neighbours(kmer.c_str(), false);
    return std::count_if(begin(nb.backward), end(nb.backward), [](auto const& lookup_result) {
        return lookup_result.kmer_id != sshash::constants::invalid_uint64;
    });
}

void DBGSSHash::call_kmers(
        const std::function<void(node_index, const std::string&)>& callback) const {
    for (size_t node_idx = 1; node_idx <= num_nodes(); ++node_idx) {
        callback(node_idx, get_node_sequence(node_idx));
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    return num_nodes() ? dict_.lookup(kmer.data(), false) + 1 : npos;
}

std::pair<DBGSSHash::node_index, bool>
DBGSSHash::kmer_to_node_with_rc(std::string_view kmer) const {
    if (!num_nodes())
        return std::make_pair(npos, false);

    auto res = dict_.lookup_advanced(kmer.data(), true);
    return std::make_pair(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    std::string str_kmer(k_, ' ');
    uint64_t ssh_idx = graph_index_to_sshash(node);
    dict_.access(ssh_idx, str_kmer.data());
    return str_kmer;
}

void DBGSSHash::serialize(std::ostream& out) const {
    // TODO
    throw std::runtime_error("serialize to stream not implemented");
}

void DBGSSHash::serialize(const std::string& filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);

    // TODO: fix this in the essentials library. for some reason, it's saver takes a non-const ref
    essentials::save(const_cast<sshash::dictionary<kmer_t>&>(dict_),
                     suffixed_filename.c_str());
}

bool DBGSSHash::load(std::istream& in) {
    throw std::runtime_error("load from stream not implemented");
    return false;
}

bool DBGSSHash::load(const std::string& filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    uint64_t num_bytes_read = essentials::load(dict_, suffixed_filename.c_str());
    if (common::get_verbose()) {
        std::cerr << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict_.size() << " [bits/kmer])"
                  << std::endl;
        dict_.print_info();
    }
    k_ = dict_.k();
    return true;
}

bool DBGSSHash::operator==(const DeBruijnGraph& other) const {
    throw std::runtime_error("operator== not implemented");
    return false;
}

const std::string DBGSSHash::alphabet_ = "ACGT";
const std::string& DBGSSHash::alphabet() const {
    return alphabet_;
}

uint64_t DBGSSHash::num_nodes() const {
    return dict_.size();
}
} // namespace graph
} // namespace mtg
