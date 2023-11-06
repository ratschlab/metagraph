#include "dbg_sshash.hpp"

namespace mtg {
namespace graph {

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)> &on_insertion) {
}

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
}

void DBGSSHash ::map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate) const {
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    return 0;
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    return 0;
}

void DBGSSHash ::adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
}

void DBGSSHash ::adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
}

void DBGSSHash ::call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const {
}

void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const {
}

size_t DBGSSHash::outdegree(node_index node) const {
    return 0;
}

bool DBGSSHash::has_single_outgoing(node_index node) const {
    return false;
}

bool DBGSSHash::has_multiple_outgoing(node_index node) const {
    return false;
}

size_t DBGSSHash::indegree(node_index node) const {
    return 0;
}

bool DBGSSHash::has_no_incoming(node_index node) const {
    return true;
}

bool DBGSSHash::has_single_incoming(node_index node) const {
    return false;
}

void DBGSSHash ::call_kmers(
        const std::function<void(node_index, const std::string &)> &callback) const {
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    return 0;
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    return "";
}

void DBGSSHash::serialize(std::ostream &out) const {
}

void DBGSSHash::serialize(const std::string &filename) const {
}

bool DBGSSHash::load(std::istream &in) {
    return false;
}

bool DBGSSHash::load(const std::string &filename) {
    return false;
}

bool DBGSSHash::operator==(const DeBruijnGraph &other) const {
    return false;
}

const std::string DBGSSHash::alphabet_ = "";
const std::string &DBGSSHash::alphabet() const {
    return alphabet_;
}

} // namespace graph
} // namespace mtg
