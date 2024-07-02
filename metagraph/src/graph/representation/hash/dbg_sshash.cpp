#include "dbg_sshash.hpp"

#include <type_traits>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/logger.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {

#if _PROTEIN_GRAPH
static constexpr uint16_t bits_per_char = sshash::aa_uint_kmer_t<uint64_t>::bits_per_char;
#else
static constexpr uint16_t bits_per_char = sshash::dna_uint_kmer_t<uint64_t>::bits_per_char;
#endif

template <typename T>
struct template_parameter;

template <template <typename ...> class C, typename T>
struct template_parameter<C<T>> { using type = T; };

template <typename T>
using get_kmer_t = typename template_parameter<std::decay_t<T>>::type;

const std::string DBGSSHash::alphabet_ = kmer::KmerExtractor2Bit().alphabet;

constexpr DeBruijnGraph::node_index sshash_to_graph_index(uint64_t idx) {
    return idx + 1;
}
constexpr uint64_t graph_index_to_sshash(DeBruijnGraph::node_index idx) {
    return idx - 1;
}

size_t DBGSSHash::dict_size() const {
    return std::visit([](const auto &d) { return d.size(); }, dict_);
}

DBGSSHash::DBGSSHash(size_t k, Mode mode) : k_(k), num_nodes_(0), mode_(mode) {
    size_t odd_k = (k_ | 1);

    if (odd_k * bits_per_char <= 64) {
        dict_ = sshash::dictionary<kmer_t<KmerInt64>>();
    } else if (odd_k * bits_per_char <= 128) {
        dict_ = sshash::dictionary<kmer_t<KmerInt128>>();
    } else if (odd_k * bits_per_char <= 256) {
        dict_ = sshash::dictionary<kmer_t<KmerInt256>>();
    } else {
        common::logger->error("k too big: {} > {}", odd_k, 256 / bits_per_char);
        throw std::length_error("k fail");
    }
}

DBGSSHash::DBGSSHash(const std::string &input_filename, size_t k, Mode mode, size_t num_chars)
      : DBGSSHash(k, mode) {
    if (k <= 1)
        throw std::domain_error("k must be at least 2");

    sshash::build_configuration build_config;
    build_config.k = k;

    // use a value of m recommended here:
    // https://twitter.com/giulio_pibiri/status/1803417277213114501
    // otherwise, use the next odd value greater than or equal to floor(k / 2)
    build_config.m = std::min({
        (num_chars > 0
            ? uint64_t(ceil(log(static_cast<double>(num_chars)) / log(static_cast<double>(alphabet_.size()))) + 1)
            : uint64_t(k / 2)) | 1,
        std::visit([](const auto &dict) -> uint64_t {
            return get_kmer_t<decltype(dict)>::max_m;
        }, dict_),
        uint64_t(k - 2)
    }) | 1;

    build_config.verbose = common::get_verbose();
    build_config.num_threads = get_num_threads();

    // silence sshash construction messages when not verbose
    std::ios orig_state(nullptr);
    orig_state.copyfmt(std::cout);
    if (!common::get_verbose())
        std::cout.setstate(std::ios_base::failbit);
    std::visit([&](auto &d) { d.build(input_filename, build_config); }, dict_);
    if (!common::get_verbose())
        std::cout.copyfmt(orig_state);

    num_nodes_ = dict_size();
}

std::string DBGSSHash::file_extension() const {
    return kExtension;
}

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)>& on_insertion) {
    throw std::logic_error("adding sequences not supported");
}

template <bool with_rc>
void DBGSSHash::map_to_nodes_with_rc(std::string_view sequence,
                                     const std::function<void(node_index, bool)>& callback,
                                     const std::function<bool()>& terminate) const {
    if (terminate() || sequence.size() < k_)
        return;

    if (!num_nodes()) {
        for (size_t i = 0; i < sequence.size() - k_ + 1 && !terminate(); ++i) {
            callback(npos, false);
        }
        return;
    }

    std::visit([&](const auto &dict) {
        using kmer_t = get_kmer_t<decltype(dict)>;
        kmer_t uint_kmer = sshash::util::string_to_uint_kmer<kmer_t>(sequence.data(), k_ - 1);
        uint_kmer.pad_char();
        for (size_t i = k_ - 1; i < sequence.size() && !terminate(); ++i) {
            uint_kmer.drop_char();
            uint_kmer.kth_char_or(k_ - 1, kmer_t::char_to_uint(sequence[i]));
            auto res = dict.lookup_advanced_uint(uint_kmer, with_rc);
            callback(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
        }
    }, dict_);
}
template
void DBGSSHash::map_to_nodes_with_rc<true>(std::string_view,
                                           const std::function<void(node_index, bool)>&,
                                           const std::function<bool()>&) const;
template
void DBGSSHash::map_to_nodes_with_rc<false>(std::string_view,
                                            const std::function<void(node_index, bool)>&,
                                            const std::function<bool()>&) const;

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)>& callback,
                             const std::function<bool()>& terminate) const {
    if (mode_ != CANONICAL) {
        map_to_nodes_sequentially(sequence, callback, terminate);
    } else {
        map_to_nodes_with_rc<true>(
                sequence, [&](node_index node, bool) { callback(node); }, terminate);
    }
}

DBGSSHash::node_index DBGSSHash::reverse_complement(node_index node) const {
    if (node == npos)
        return npos;

    assert(node > 0 && node <= num_nodes());
    if (node > dict_size())
        return node - dict_size();

    if (k_ % 2 == 1)
        return node + dict_size();

    std::string str_kmer(k_, ' ');
    uint64_t ssh_idx = graph_index_to_sshash(node);
    std::visit([&](const auto &d) { d.access(ssh_idx, str_kmer.data()); }, dict_);
    std::string rc_str_kmer(str_kmer);
    ::reverse_complement(rc_str_kmer.begin(), rc_str_kmer.end());
    return str_kmer != rc_str_kmer ? node + dict_size() : node;
}

uint64_t DBGSSHash::num_nodes() const {
    return mode_ != CANONICAL ? dict_size() : dict_size() * 2;
}

void DBGSSHash::map_to_nodes_sequentially(std::string_view sequence,
                                          const std::function<void(node_index)>& callback,
                                          const std::function<bool()>& terminate) const {
    if (mode_ == CANONICAL) {
        map_to_nodes_with_rc<true>(sequence, [&](node_index n, bool orientation) {
            callback(orientation ? reverse_complement(n) : n);
        }, terminate);
    } else {
        map_to_nodes_with_rc<false>(sequence, [&](node_index node, bool) {
            callback(node);
        }, terminate);
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    assert(node > 0 && node <= num_nodes());
    // TODO: if a node is in the middle of a unitig, then we only need to check the next node index
    return kmer_to_node(get_node_sequence(node).substr(1) + next_char);
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    assert(node > 0 && node <= num_nodes());
    // TODO: if a node is in the middle of a unitig, then we only need to check the previous node index
    std::string string_kmer = prev_char + get_node_sequence(node);
    string_kmer.pop_back();
    return kmer_to_node(string_kmer);
}

template <bool with_rc>
void DBGSSHash::call_outgoing_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    std::visit([&](const auto &dict) {
        using kmer_t = get_kmer_t<decltype(dict)>;
        auto nb = dict.kmer_forward_neighbours(kmer.c_str(), with_rc);
        for (size_t i = 0; i < nb.forward.size(); i++) {
            if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                callback(sshash_to_graph_index(nb.forward[i].kmer_id),
                         kmer_t::uint64_to_char(i),
                         nb.forward[i].kmer_orientation);
            }
        }
    }, dict_);
}
template
void DBGSSHash::call_outgoing_kmers_with_rc<true>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;
template
void DBGSSHash::call_outgoing_kmers_with_rc<false>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;

template <bool with_rc>
void DBGSSHash::call_incoming_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    std::visit([&](const auto &dict) {
        using kmer_t = get_kmer_t<decltype(dict)>;
        auto nb = dict.kmer_backward_neighbours(kmer.c_str(), with_rc);
        for (size_t i = 0; i < nb.backward.size(); i++) {
            if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                callback(sshash_to_graph_index(nb.backward[i].kmer_id),
                         kmer_t::uint64_to_char(i),
                         nb.backward[i].kmer_orientation);
            }
        }
    }, dict_);
}
template
void DBGSSHash::call_incoming_kmers_with_rc<true>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;
template
void DBGSSHash::call_incoming_kmers_with_rc<false>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;

void DBGSSHash::call_outgoing_kmers(node_index node,
                                    const OutgoingEdgeCallback& callback) const {
    if (mode_ == CANONICAL) {
        call_outgoing_kmers_with_rc<true>(node, [&](node_index next, char c, bool orientation) {
            callback(orientation ? reverse_complement(next) : next, c);
        });
    } else {
        call_outgoing_kmers_with_rc<false>(node, [&](node_index next, char c, bool) {
            callback(next, c);
        });
    }
}

void DBGSSHash::call_incoming_kmers(node_index node,
                                    const IncomingEdgeCallback& callback) const {
    if (mode_ == CANONICAL) {
        call_incoming_kmers_with_rc<true>(node, [&](node_index prev, char c, bool orientation) {
            callback(orientation ? reverse_complement(prev) : prev, c);
        });
    } else {
        call_incoming_kmers_with_rc<false>(node, [&](node_index prev, char c, bool) {
            callback(prev, c);
        });
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    size_t res = 0;
    adjacent_outgoing_nodes(node, [&](node_index) { ++res; });
    return res;
}

size_t DBGSSHash::indegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    size_t res = 0;
    adjacent_incoming_nodes(node, [&](node_index) { ++res; });
    return res;
}

void DBGSSHash::call_nodes(
        const std::function<void(node_index)>& callback,
        const std::function<bool()> &terminate) const {
    for (size_t node_idx = 1; !terminate() && node_idx <= dict_size(); ++node_idx) {
        callback(node_idx);
    }

    if (mode_ == CANONICAL) {
        for (size_t node_idx = 1; !terminate() && node_idx <= dict_size(); ++node_idx) {
            size_t rc_node_idx = reverse_complement(node_idx);
            if (rc_node_idx != node_idx)
                callback(rc_node_idx);
        }
    }
}

template <bool with_rc>
std::pair<DBGSSHash::node_index, bool>
DBGSSHash::kmer_to_node_with_rc(std::string_view kmer) const {
    if (!num_nodes())
        return std::make_pair(npos, false);

    return std::visit([&](const auto &d) {
        auto res = d.lookup_advanced(kmer.data(), with_rc);
        return std::make_pair(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
    }, dict_);
}
template
std::pair<DBGSSHash::node_index, bool>
DBGSSHash::kmer_to_node_with_rc<true>(std::string_view) const;
template
std::pair<DBGSSHash::node_index, bool>
DBGSSHash::kmer_to_node_with_rc<false>(std::string_view) const;

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    if (mode_ == CANONICAL) {
        auto res = kmer_to_node_with_rc<true>(kmer);
        return res.second ? reverse_complement(res.first) : res.first;
    } else {
        return kmer_to_node_with_rc<false>(kmer).first;
    }
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    std::string str_kmer(k_, ' ');
    node_index node_canonical = node > dict_size() ? node - dict_size() : node;
    uint64_t ssh_idx = graph_index_to_sshash(node_canonical);
    std::visit([&](const auto &d) { d.access(ssh_idx, str_kmer.data()); }, dict_);

    if (node > node_canonical)
        ::reverse_complement(str_kmer.begin(), str_kmer.end());

    return str_kmer;
}

void DBGSSHash::serialize(std::ostream& out) const {
    essentials::generic_saver saver(out);

    saver.visit(num_nodes_);
    saver.visit(k_);
    saver.visit(mode_);

    if (num_nodes())
        std::visit([&](const auto &d) { saver.visit(d); }, dict_);
}

void DBGSSHash::serialize(const std::string& filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    std::ofstream fout(suffixed_filename, std::ios::binary);
    serialize(fout);
}

bool DBGSSHash::load(std::istream &in) {
    essentials::generic_loader loader(in);
    size_t num_nodes;
    size_t k;
    Mode mode;
    loader.visit(num_nodes);
    loader.visit(k);
    loader.visit(mode);

    *this = DBGSSHash(k, mode);
    num_nodes_ = num_nodes;

    if (num_nodes_)
        std::visit([&](auto &d) { d.visit(loader); }, dict_);

    return true;
}

bool DBGSSHash::load(const std::string& filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    std::ifstream fin(suffixed_filename, std::ios::binary);
    return load(fin);
}

} // namespace graph
} // namespace mtg
