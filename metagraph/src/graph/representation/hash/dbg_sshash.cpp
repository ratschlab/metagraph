#include "dbg_sshash.hpp"

#include <sdsl/uint256_t.hpp>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/logger.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {

const std::string DBGSSHash::alphabet_ = kmer::KmerExtractor2Bit().alphabet;

using Kmer64 = uint64_t;
using Kmer128 = __uint128_t;
using Kmer256 = sdsl::uint256_t;

template <typename KmerInt>
class Dictionary : public IDictionary {
  public:
#if _PROTEIN_GRAPH
    using kmer_t = sshash::aa_uint_kmer_t<KmerInt>;
#else
    using kmer_t = sshash::dna_uint_kmer_t<KmerInt>;
#endif

    const sshash::dictionary<kmer_t>& data() const { return dict_; }
    size_t size() const { return dict_.size(); }

    void visit(essentials::loader &loader) { dict_.visit(loader); }
    void visit(essentials::saver &saver) const { dict_.visit(saver); }

    void access(uint64_t kmer_id, char* string_kmer) const { dict_.access(kmer_id, string_kmer); }
    void build(std::string const& input_filename, sshash::build_configuration const& build_config) {
        dict_.build(input_filename, build_config);
    }

    uint64_t lookup(char const* string_kmer, bool check_reverse_complement = true) const {
        return dict_.lookup(string_kmer, check_reverse_complement);
    }

    uint64_t lookup_uint(kmer_t uint_kmer, bool check_reverse_complement = true) const {
        return dict_.lookup_uint(uint_kmer, check_reverse_complement);
    }

    sshash::lookup_result lookup_advanced(
            char const* string_kmer,
            bool check_reverse_complement = true) const {
        return dict_.lookup_advanced(string_kmer, check_reverse_complement);
    }

    sshash::lookup_result lookup_advanced_uint(
            kmer_t uint_kmer,
            bool check_reverse_complement = true) const {
        return dict_.lookup_advanced_uint(uint_kmer, check_reverse_complement);
    }

    sshash::neighbourhood<kmer_t> kmer_forward_neighbours(
            char const* string_kmer,
            bool check_reverse_complement = true) const {
        return dict_.kmer_forward_neighbours(string_kmer, check_reverse_complement);
    }

    sshash::neighbourhood<kmer_t> kmer_backward_neighbours(
            char const* string_kmer,
            bool check_reverse_complement = true) const {
        return dict_.kmer_backward_neighbours(string_kmer, check_reverse_complement);
    }

  private:
    // TODO: this is a hack since sshash::dictionary<kmer_t>::visit is always non-const
    mutable sshash::dictionary<kmer_t> dict_;
};

constexpr DeBruijnGraph::node_index sshash_to_graph_index(uint64_t idx) {
    return idx + 1;
}
constexpr uint64_t graph_index_to_sshash(DeBruijnGraph::node_index idx) {
    return idx - 1;
}

DBGSSHash::DBGSSHash(size_t k, Mode mode) : k_(k), num_nodes_(0), mode_(mode) {
    size_t odd_k = (k_ | 1);

    if (odd_k * IDictionary::bits_per_char <= 64) {
        dict_ = std::make_unique<Dictionary<Kmer64>>();
    } else if (odd_k * IDictionary::bits_per_char <= 128) {
        dict_ = std::make_unique<Dictionary<Kmer128>>();
    } else if (odd_k * IDictionary::bits_per_char <= 256) {
        dict_ = std::make_unique<Dictionary<Kmer256>>();
    } else {
        common::logger->error("k too big: {} > {}", odd_k, 256 / IDictionary::bits_per_char);
        throw std::length_error("k fail");
    }
}

DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k, Mode mode)
      : DBGSSHash(k, mode) {
    sshash::build_configuration build_config;
    build_config.k = k;
    // quick fix for value of m... k/2 but odd
    build_config.m = std::min(uint64_t((k_ + 1) / 2), sshash::constants::max_m);
    build_config.verbose = common::get_verbose();
    build_config.num_threads = get_num_threads();
    if (build_config.m % 2 == 0)
        build_config.m++;

    // silence sshash construction messages when not verbose
    if (!common::get_verbose())
        std::cout.setstate(std::ios_base::failbit);

    dict_->build(input_filename, build_config);
    if (!common::get_verbose())
        std::cout.clear();

    num_nodes_ = dict_->size();
}

std::string DBGSSHash::file_extension() const {
    return kExtension;
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

DBGSSHash::node_index DBGSSHash::reverse_complement(node_index node) const {
    if (node == npos)
        return npos;

    if (node > dict_->size())
        return node - dict_->size();

    if (k_ % 2 == 1)
        return node + dict_->size();

    std::string str_kmer(k_, ' ');
    uint64_t ssh_idx = graph_index_to_sshash(node);
    dict_->access(ssh_idx, str_kmer.data());
    std::string rc_str_kmer(str_kmer);
    ::reverse_complement(rc_str_kmer.begin(), rc_str_kmer.end());
    return str_kmer != rc_str_kmer ? node + dict_->size() : node;
}

uint64_t DBGSSHash::num_nodes() const {
    return mode_ != CANONICAL ? dict_->size() : dict_->size() * 2;
}

void DBGSSHash::map_to_nodes_sequentially(std::string_view sequence,
                                          const std::function<void(node_index)>& callback,
                                          const std::function<bool()>& terminate) const {
    if (terminate() || sequence.size() < k_)
        return;

    if (!num_nodes()) {
        for (size_t i = 0; i < sequence.size() - k_ + 1 && !terminate(); ++i) {
            callback(npos);
        }
        return;
    }

    if (mode_ == CANONICAL) {
        map_to_nodes_with_rc(sequence, [&](node_index n, bool orientation) {
            callback(orientation ? reverse_complement(n) : n);
        }, terminate);
        return;
    }

    auto stream_map = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        kmer_t uint_kmer = sshash::util::string_to_uint_kmer<kmer_t>(sequence.data(), k_ - 1);
        uint_kmer.pad_char();
        for (size_t i = k_ - 1; i < sequence.size() && !terminate(); ++i) {
            uint_kmer.drop_char();
            uint_kmer.kth_char_or(k_ - 1, kmer_t::char_to_uint(sequence[i]));
            callback(sshash_to_graph_index(dict.lookup_uint(uint_kmer, false)));
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        stream_map(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        stream_map(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        stream_map(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

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

    auto stream_map = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        kmer_t uint_kmer = sshash::util::string_to_uint_kmer<kmer_t>(sequence.data(), k_ - 1);
        uint_kmer.pad_char();
        for (size_t i = k_ - 1; i < sequence.size() && !terminate(); ++i) {
            uint_kmer.drop_char();
            uint_kmer.kth_char_or(k_ - 1, kmer_t::char_to_uint(sequence[i]));
            auto res = dict.lookup_advanced_uint(uint_kmer, true);
            callback(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        stream_map(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        stream_map(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        stream_map(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    // TODO: if a node is in the middle of a unitig, then we only need to check the next node index
    return kmer_to_node(get_node_sequence(node).substr(1) + next_char);
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    // TODO: if a node is in the middle of a unitig, then we only need to check the previous node index
    std::string string_kmer = std::string(1, prev_char) + get_node_sequence(node);
    string_kmer.pop_back();
    return kmer_to_node(string_kmer);
}

void DBGSSHash::adjacent_outgoing_nodes(node_index node,
                                        const std::function<void(node_index)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

void DBGSSHash::adjacent_incoming_nodes(node_index node,
                                        const std::function<void(node_index)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGSSHash::call_outgoing_kmers(node_index node,
                                    const OutgoingEdgeCallback& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    auto call_impl = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        sshash::neighbourhood<kmer_t> nb = dict.kmer_forward_neighbours(kmer.c_str(), mode_ == CANONICAL);
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                node_index next = sshash_to_graph_index(nb.forward[i].kmer_id);
                callback(nb.forward[i].kmer_orientation ? reverse_complement(next) : next, kmer_t::alphabet[i]);
            }
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        call_impl(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        call_impl(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        call_impl(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

void DBGSSHash::call_incoming_kmers(node_index node,
                                    const IncomingEdgeCallback& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    auto call_impl = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        sshash::neighbourhood<kmer_t> nb = dict.kmer_backward_neighbours(kmer.c_str(), mode_ == CANONICAL);
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                node_index prev = sshash_to_graph_index(nb.backward[i].kmer_id);
                callback(nb.backward[i].kmer_orientation ? reverse_complement(prev) : prev, kmer_t::alphabet[i]);
            }
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        call_impl(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        call_impl(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        call_impl(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

void DBGSSHash::call_outgoing_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    auto call_impl = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        sshash::neighbourhood<kmer_t> nb = dict.kmer_forward_neighbours(kmer.c_str(), true);
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                callback(sshash_to_graph_index(nb.forward[i].kmer_id), kmer_t::alphabet[i],
                         nb.forward[i].kmer_orientation);
            }
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        call_impl(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        call_impl(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        call_impl(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

void DBGSSHash::call_incoming_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(node > 0 && node <= num_nodes());
    std::string kmer = get_node_sequence(node);
    auto call_impl = [&](const auto &dict) {
        using kmer_t = typename std::remove_reference<decltype(dict)>::type::kmer_t;
        sshash::neighbourhood<kmer_t> nb = dict.kmer_backward_neighbours(kmer.c_str(), true);
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                callback(sshash_to_graph_index(nb.backward[i].kmer_id), kmer_t::alphabet[i],
                         nb.backward[i].kmer_orientation);
            }
        }
    };

    if (auto dict64 = dynamic_cast<Dictionary<Kmer64>*>(dict_.get())) {
        call_impl(*dict64);
    } else if (auto dict128 = dynamic_cast<Dictionary<Kmer128>*>(dict_.get())) {
        call_impl(*dict128);
    } else if (auto dict256 = dynamic_cast<Dictionary<Kmer256>*>(dict_.get())) {
        call_impl(*dict256);
    } else {
        common::logger->error("");
        throw std::runtime_error("k fail");
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    size_t res = 0;
    adjacent_outgoing_nodes(node, [&](node_index) { ++res; });
    return res;
}

size_t DBGSSHash::indegree(node_index node) const {
    size_t res = 0;
    adjacent_incoming_nodes(node, [&](node_index) { ++res; });
    return res;
}

void DBGSSHash::call_nodes(
        const std::function<void(node_index)>& callback,
        const std::function<bool()> &terminate) const {
    for (size_t node_idx = 1; !terminate() && node_idx <= dict_->size(); ++node_idx) {
        callback(node_idx);
        if (mode_ == CANONICAL && !terminate()) {
            size_t rc_node_idx = reverse_complement(node_idx);
            if (rc_node_idx != node_idx)
                callback(rc_node_idx);
        }
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    if (!num_nodes())
        return npos;

    auto res = dict_->lookup_advanced(kmer.data(), mode_ == CANONICAL);
    node_index node = sshash_to_graph_index(res.kmer_id);
    return res.kmer_orientation ? reverse_complement(node) : node;
}

std::pair<DBGSSHash::node_index, bool>
DBGSSHash::kmer_to_node_with_rc(std::string_view kmer) const {
    if (!num_nodes())
        return std::make_pair(npos, false);

    auto res = dict_->lookup_advanced(kmer.data(), true);
    return std::make_pair(sshash_to_graph_index(res.kmer_id), res.kmer_orientation);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    std::string str_kmer(k_, ' ');
    node_index node_canonical = node > dict_->size() ? node - dict_->size() : node;
    uint64_t ssh_idx = graph_index_to_sshash(node_canonical);
    dict_->access(ssh_idx, str_kmer.data());

    if (node > node_canonical)
        ::reverse_complement(str_kmer.begin(), str_kmer.end());

    return str_kmer;
}

void DBGSSHash::serialize(std::ostream& out) const {
    // TODO
    throw std::runtime_error("serialize to stream not implemented");
}

void DBGSSHash::serialize(const std::string& filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    essentials::saver saver(suffixed_filename.c_str());

    // TODO: fix this in the essentials library. for some reason, its saver takes a non-const ref
    saver.visit(const_cast<size_t&>(num_nodes_));
    saver.visit(const_cast<size_t&>(k_));
    saver.visit(const_cast<Mode&>(mode_));

    if (num_nodes())
        dict_->visit(saver);
}

bool DBGSSHash::load(std::istream& in) {
    throw std::runtime_error("load from stream not implemented");
    return false;
}

bool DBGSSHash::load(const std::string& filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    essentials::loader loader(suffixed_filename.c_str());

    size_t num_nodes;
    size_t k;
    Mode mode;
    loader.visit(num_nodes);
    loader.visit(k);
    loader.visit(mode);

    *this = DBGSSHash(k, mode);
    num_nodes_ = num_nodes;

    if (num_nodes_)
        dict_->visit(loader);

    return true;
}

bool DBGSSHash::operator==(const DeBruijnGraph& other) const {
    try {
        other.call_kmers([&](node_index, const std::string &kmer) {
            if (!find(kmer))
                throw std::bad_function_call();
        });

        call_kmers([&](node_index, const std::string &kmer) {
            if (!other.find(kmer))
                throw std::bad_function_call();
        });
    } catch (const std::bad_function_call&) {
        return false;
    }

    return true;
}

} // namespace graph
} // namespace mtg
