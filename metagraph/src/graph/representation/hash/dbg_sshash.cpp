#include "dbg_sshash.hpp"
#include <dictionary.hpp>


namespace mtg {
namespace graph {
DBGSSHash::~DBGSSHash() {}

DBGSSHash::DBGSSHash(size_t k):k_(k) {
    dict_ = std::make_unique<sshash::dictionary>();
}

DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k):k_(k){
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
    DeBruijnGraph::Mode DBGSSHash::get_mode() const  { return BASIC; }

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
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
        callback(kmer_to_node(sequence.substr(i, k_)));
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node); 
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0]);
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
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(&kmer[0]);
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

    auto prefix = get_node_sequence(node).substr(1);

    for (char c : alphabet_) {
        auto next = kmer_to_node(prefix + c);
        if (next != npos)
            callback(next, c);
    }
}


void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string suffix = get_node_sequence(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto prev = kmer_to_node(c + suffix);
        if (prev != npos)
            callback(prev, c);
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0]);
    size_t out_deg = bool(nb.forward_A.kmer_id + 1) // change to loop?
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
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(&kmer[0]);
    size_t in_deg = bool(nb.backward_A.kmer_id + 1) // change to loop?
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
    for (size_t node_idx = 1; node_idx <= num_nodes(); ++node_idx) {
        callback(node_idx, get_node_sequence(node_idx));
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    uint64_t ssh_idx = dict_->lookup(kmer.begin(), false);
    return ssh_idx + 1;
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
    
    essentials::logger("saving data structure to disk...");
    essentials::save(*dict_, suffixed_filename.c_str());
    essentials::logger("DONE");
}

bool DBGSSHash::load(std::istream &in) {
    throw std::runtime_error("load from stream not implemented");
    return false;
}

bool DBGSSHash::load(const std::string &filename) {

    uint64_t num_bytes_read = essentials::load(*dict_, filename.c_str());
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
