#include "dbg_sshash.hpp"
#include <dictionary.hpp>

#include <query/streaming_query_regular_parsing.hpp>


namespace mtg {
namespace graph {
DBGSSHash::~DBGSSHash() {}

DBGSSHash::DBGSSHash(size_t k):k_(k) {
    dict_ = std::make_unique<sshash::dictionary>();
}

DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k, Mode mode):k_(k), mode_(mode) {
    sshash::build_configuration build_config;
    build_config.k = k;//
    // quick fix for value of m... k/2 but odd
    build_config.m = (k_+1)/2;
    if(build_config.m % 2 == 0) build_config.m++;
    dict_ = std::make_unique<sshash::dictionary>();
    dict_->build(input_filename, build_config);
}
    std::string DBGSSHash::file_extension() const  { return kExtension; }
    size_t DBGSSHash::get_k() const  { return k_; }
    DeBruijnGraph::Mode DBGSSHash::get_mode() const  { return mode_; }

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)> &on_insertion) {
    // TODO: throw exception? :)
        throw std::runtime_error("adding sequences not implemented");

}

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    map_to_nodes_sequentially(sequence, callback, terminate);
}

void DBGSSHash ::map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate) const {
    if (terminate() || sequence.size() < k_)
        return;

    auto uint_kmer = sshash::util::string_to_uint_kmer(sequence.data(), k_ - 1) << 2;
    for (size_t i = k_ - 1; i < sequence.size() && !terminate(); ++i) {
        uint_kmer = (uint_kmer >> 2) + (sshash::util::char_to_uint(sequence[i]) << (2 * (k_ - 1)));
        callback(dict_->lookup_uint(uint_kmer, false) + 1);
    }
}

void DBGSSHash ::map_to_nodes_with_rc(std::string_view sequence,
                                      const std::function<void(node_index, bool)> &callback,
                                      const std::function<bool()> &terminate) const {
    sshash::streaming_query_regular_parsing streamer(dict_.get());
    streamer.start();
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
        const char *kmer = sequence.data() + i;
        auto res = streamer.lookup_advanced(kmer);
        callback(res.kmer_id + 1, res.kmer_orientation);
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0], false);
    uint64_t ssh_idx = -1;
    switch (next_char) {
    case 'A':
        ssh_idx = nb.forward_A.kmer_id;
        break;
    case 'C':
        ssh_idx = nb.forward_C.kmer_id;
        break;
    case 'G':
        ssh_idx = nb.forward_G.kmer_id;
        break;
    case 'T':
        ssh_idx = nb.forward_T.kmer_id;
        break;
    default:
        break;
    }
    return ssh_idx + 1;
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(&kmer[0], false);
    uint64_t ssh_idx = -1;
    switch (prev_char) {
    case 'A':
        ssh_idx = nb.backward_A.kmer_id;
        break;
    case 'C':
        ssh_idx = nb.backward_C.kmer_id;
        break;
    case 'G':
        ssh_idx = nb.backward_G.kmer_id;
        break;
    case 'T':
        ssh_idx = nb.backward_T.kmer_id;
        break;
    default:
        break;
    }
    return ssh_idx + 1;
}

void DBGSSHash ::adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());
    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

void DBGSSHash ::adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());
    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGSSHash ::call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(kmer.c_str(), false);
    if (nb.forward_A.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_A.kmer_id + 1, 'A');

    if (nb.forward_C.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_C.kmer_id + 1, 'C');

    if (nb.forward_G.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_G.kmer_id + 1, 'G');

    if (nb.forward_T.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_T.kmer_id + 1, 'T');
}


void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(kmer.c_str(), false);
    if (nb.backward_A.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_A.kmer_id + 1, 'A');

    if (nb.backward_C.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_C.kmer_id + 1, 'C');

    if (nb.backward_G.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_G.kmer_id + 1, 'G');

    if (nb.backward_T.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_T.kmer_id + 1, 'T');
}

void DBGSSHash ::call_outgoing_kmers_with_rc(node_index node,
                                     const std::function<void(node_index, char, bool)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(kmer.c_str(), true);
    if (nb.forward_A.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_A.kmer_id + 1, 'A', nb.forward_A.kmer_orientation);

    if (nb.forward_C.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_C.kmer_id + 1, 'C', nb.forward_C.kmer_orientation);

    if (nb.forward_G.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_G.kmer_id + 1, 'G', nb.forward_G.kmer_orientation);

    if (nb.forward_T.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.forward_T.kmer_id + 1, 'T', nb.forward_T.kmer_orientation);
}


void DBGSSHash ::call_incoming_kmers_with_rc(node_index node,
                                     const std::function<void(node_index, char, bool)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(kmer.c_str(), true);
    if (nb.backward_A.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_A.kmer_id + 1, 'A', nb.backward_A.kmer_orientation);

    if (nb.backward_C.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_C.kmer_id + 1, 'C', nb.backward_C.kmer_orientation);

    if (nb.backward_G.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_G.kmer_id + 1, 'G', nb.backward_G.kmer_orientation);

    if (nb.backward_T.kmer_id != sshash::constants::invalid_uint64)
        callback(nb.backward_T.kmer_id + 1, 'T', nb.backward_T.kmer_orientation);
}

size_t DBGSSHash::outdegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0], false);
    size_t out_deg = (nb.forward_A.kmer_id != sshash::constants::invalid_uint64) // change to loop?
                    + (nb.forward_C.kmer_id != sshash::constants::invalid_uint64)
                    + (nb.forward_G.kmer_id != sshash::constants::invalid_uint64)
                    + (nb.forward_T.kmer_id != sshash::constants::invalid_uint64);
    return out_deg;
}

bool DBGSSHash::has_single_outgoing(node_index node) const {
    return DBGSSHash::outdegree(node) == 1;
}

bool DBGSSHash::has_multiple_outgoing(node_index node) const {
    return DBGSSHash::outdegree(node) > 1;
}

size_t DBGSSHash::indegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(kmer.c_str(), false);
    size_t in_deg = (nb.backward_A.kmer_id != sshash::constants::invalid_uint64) // change to loop?
                    + (nb.backward_C.kmer_id != sshash::constants::invalid_uint64)
                    + (nb.backward_G.kmer_id != sshash::constants::invalid_uint64)
                    + (nb.backward_T.kmer_id != sshash::constants::invalid_uint64);
    return in_deg;
}

bool DBGSSHash::has_no_incoming(node_index node) const {
    return DBGSSHash::indegree(node) == 0;
}

bool DBGSSHash::has_single_incoming(node_index node) const {
    return DBGSSHash::indegree(node) == 1;
}

void DBGSSHash::call_kmers(
        const std::function<void(node_index, const std::string &)> &callback) const {
    for (size_t node_idx = 1; node_idx <= num_nodes(); ++node_idx) {
        callback(node_idx, get_node_sequence(node_idx));
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    return num_nodes() ? dict_->lookup(kmer.begin(), false) + 1 : npos;
}

std::pair<DBGSSHash::node_index, bool> DBGSSHash::kmer_to_node_with_rc(std::string_view kmer) const {
    if (!num_nodes())
        return std::make_pair(npos, false);

    auto res = dict_->lookup_advanced(kmer.begin(), true);
    return std::make_pair(res.kmer_id + 1, res.kmer_orientation);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    std::string str_kmer = "";
    str_kmer.append(k_, ' ');
    uint64_t ssh_idx = node - 1; // switch back to sshash idx!!!
    dict_->access(ssh_idx, &str_kmer[0]);
    return str_kmer;
}

void DBGSSHash::serialize(std::ostream &out) const {
    //TODO
    throw std::runtime_error("serialize to stream not implemented");
}

void DBGSSHash::serialize(const std::string &filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    
    common::logger->trace("saving data structure to disk...");
    essentials::save(*dict_, suffixed_filename.c_str());
    essentials::logger("DONE");
}

bool DBGSSHash::load(std::istream &in) {
    throw std::runtime_error("load from stream not implemented");
    return false;
}

bool DBGSSHash::load(const std::string &filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    uint64_t num_bytes_read = essentials::load(*dict_, suffixed_filename.c_str());
    bool verbose = true; // temp
    if (verbose) {
        std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict_->size() << " [bits/kmer])"
                  << std::endl;
        dict_->print_info();
    }
    k_ = dict_->k();
    return true;
}

bool DBGSSHash::operator==(const DeBruijnGraph &other) const {
    throw std::runtime_error("operator== not implemented");
    return false;
}

const std::string DBGSSHash::alphabet_ = "ACGT";
const std::string &DBGSSHash::alphabet() const {
    return alphabet_;
}

uint64_t DBGSSHash::num_nodes() const { return dict_->size(); }
} // namespace graph
} // namespace mtg
