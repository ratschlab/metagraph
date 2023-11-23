#include "dbg_sshash.hpp"

namespace mtg {
namespace graph {

// replace this function by just using other constructor and load?
DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k):k_(k){
    sshash::build_configuration build_config;
    build_config.k = k;
    dict_.build(input_filename, build_config);
}
void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)> &on_insertion) {
    // TODO: throw exception? :)
}

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    map_to_nodes_sequentially(sequence, callback, terminate);
}

void DBGSSHash ::map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate) const {
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
        callback(kmer_to_node(sequence.substr(i, k_)));
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_.kmer_forward_neighbours(&kmer[0]);
    switch (next_char) {
    case 'A':
        return nb.forward_A.kmer_id;
    case 'C':
        return nb.forward_C.kmer_id;
    case 'G':
        return nb.forward_G.kmer_id;
    case 'T':
        return nb.forward_T.kmer_id;
    default:
        break;
    }
    return 0;
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_.kmer_backward_neighbours(&kmer[0]);
    switch (prev_char) {
    case 'A':
        return nb.backward_A.kmer_id;
    case 'C':
        return nb.backward_C.kmer_id;
    case 'G':
        return nb.backward_G.kmer_id;
    case 'T':
        return nb.backward_T.kmer_id;
    default:
        break;
    }
    return 0;
}

void DBGSSHash ::adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node+1 > 0 && node < num_nodes());
    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
    
}

void DBGSSHash ::adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node+1 > 0 && node < num_nodes());
    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGSSHash ::call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const {
    assert(node+1 > 0 && node < num_nodes());

    auto prefix = get_node_sequence(node).substr(1);

    for (char c : alphabet_) {
        auto next = kmer_to_node(prefix + c);
        if (next+1 != npos)
            // ! next ranges from 0 to num_nodes -1 !
            callback(next, c);
    }
}


void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const {
    assert(node+1 > 0 && node < num_nodes());

    std::string suffix = get_node_sequence(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto prev = kmer_to_node(c + suffix);
        if (prev+1 != npos)
            // ! prev ranges from 0 to num_nodes -1 !
            callback(prev, c);
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_.kmer_forward_neighbours(&kmer[0]);
    size_t out_deg = bool(nb.forward_A.kmer_id + 1)
                    + bool(nb.forward_C.kmer_id + 1)
                    + bool(nb.forward_G.kmer_id + 1)
                    + bool(nb.forward_T.kmer_id + 1);
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
    sshash::neighbourhood nb = dict_.kmer_backward_neighbours(&kmer[0]);
    size_t in_deg = bool(nb.backward_A.kmer_id + 1) 
                    + bool(nb.backward_C.kmer_id + 1)
                    + bool(nb.backward_G.kmer_id + 1)
                    + bool(nb.backward_T.kmer_id + 1);
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
    for (size_t node_idx = 0; node_idx < num_nodes(); ++node_idx) {
        callback(node_idx, get_node_sequence(node_idx));
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    return dict_.lookup(kmer.begin(), true);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    std::string str_kmer = "";
    str_kmer.append(k_, ' ');
    dict_.access(node, &str_kmer[0]);
    return str_kmer;
}

void DBGSSHash::serialize(std::ostream &out) const {
    //TODO
    throw std::invalid_argument( "not implemented" );
}

void DBGSSHash::serialize(const std::string &filename) const {
    dict_.dump(filename);
}

bool DBGSSHash::load(std::istream &in) {
    //TODO
    return false;
}

bool DBGSSHash::load(const std::string &filename) {
    sshash::build_configuration build_config;
    build_config.k = k_;
    // quick fix for value of m... k/2
    build_config.m = (k_+1)/2;
    //build_config.tmp_dirname = "/home/marianna/Documents/Masterthesis/metagraph/metagraph/tests/data/sshash_sequences";
    dict_.build(filename, build_config);
    return true;
}

bool DBGSSHash::operator==(const DeBruijnGraph &other) const {
    //TODO
    return false;
}

const std::string DBGSSHash::alphabet_ = "ACGT";
const std::string &DBGSSHash::alphabet() const {
    return alphabet_;
}

} // namespace graph
} // namespace mtg
