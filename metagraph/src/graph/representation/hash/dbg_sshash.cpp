#include "dbg_sshash.hpp"

#include <algorithm>
#include <bit>
#include <mutex>
#include <type_traits>

#include <sdsl/int_vector.hpp>
#include <sdsl/memory_management.hpp>

#include <streaming_query.hpp>
#include <progress_bar.hpp>
#include <spectrum_preserving_string_set.hpp>
#include <kmer_iterator.hpp>

// sshash template implementations
#include <src/dictionary.cpp>
#include <src/builder/build.cpp>
#include <src/query.cpp>
#include <src/info.cpp>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/file_utils.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {

#if _PROTEIN_GRAPH
static constexpr uint16_t bits_per_char = sshash::aa_uint_kmer_t<uint64_t>::bits_per_char;
#else
static constexpr uint16_t bits_per_char = sshash::dna_uint_kmer_t<uint64_t>::bits_per_char;
#endif


#if _PROTEIN_GRAPH
const std::string DBGSSHash::alphabet_ = kmer::alphabets::kAlphabetProtein;
#else
// SSHash's nucleotide dictionary has a two-bit alphabet, including in a DNA5
// build. N must interrupt lookup instead of being encoded as a fifth symbol.
const std::string DBGSSHash::alphabet_ = kmer::alphabets::kAlphabetDNA;
#endif

// Only node IDs survive construction, so copies may share the order without
// retaining references to any particular graph or dictionary.
struct DBGSSHash::SeedIndexCache {
    std::once_flag initialized;
    std::once_flag suffix_initialized;
    sdsl::int_vector<> prefix_order;
    sdsl::int_vector<> suffix_order;

    template <bool prefix, class Dictionary>
    void build(const DBGSSHash& graph, const Dictionary& dictionary, sdsl::int_vector<>& order) {
        using Kmer = typename Dictionary::kmer_type;
        using Word = decltype(Kmer{}.bits);
        const size_t k = graph.get_k();
        const bool canonical = graph.get_mode() == CANONICAL;
        if (graph.get_mode() != BASIC && !canonical)
            throw std::invalid_argument("Seed lookup requires BASIC or CANONICAL SSHash");
        if constexpr (!Kmer::has_reverse_complement) {
            if (canonical)
                throw std::invalid_argument("Canonical seed lookup requires reverse complements");
        }

        std::vector<std::pair<Word, node_index>> entries;
        entries.reserve(graph.num_nodes());
        for (auto it = dictionary.begin(); it.has_next();) {
            auto [id, word] = it.next();
            // Low characters come first in the packed word. Prefix lookup
            // reverses them; suffix lookup keeps their colexicographic order.
            Kmer key = word;
            if constexpr (prefix) {
                key = Kmer{};
                for (size_t i = 0; i < k; ++i)
                    key.set(k - i - 1, word.at(i));
            }
            entries.emplace_back(key.bits, sshash_to_graph_index(id));
            if (canonical) {
                const auto reverse_node = graph.reverse_complement(sshash_to_graph_index(id));
                if (reverse_node != sshash_to_graph_index(id)) {
                    key.reverse_complement_inplace(k);
                    entries.emplace_back(key.bits, reverse_node);
                }
            }
        }
        // Pair's C++20 comparison can narrow uint256_t through its uint64_t conversion.
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            return a.first == b.first ? a.second < b.second : a.first < b.first;
        });
        const auto width = std::max<unsigned>(1, std::bit_width(graph.max_index()));
        order = sdsl::int_vector<>(entries.size(), 0, width);
        for (size_t i = 0; i < entries.size(); ++i)
            order[i] = entries[i].second;
    }

    template <class Dictionary>
    static auto word_at(const DBGSSHash& graph, const Dictionary& dictionary, node_index node) {
        const size_t base_nodes = dictionary.num_kmers();
        const bool reverse = node > base_nodes;
        const auto id = graph_index_to_sshash(node) - (reverse ? base_nodes : 0);
        auto word = dictionary.at_kmer_id(id).next().second;
        if (reverse) word.reverse_complement_inplace(graph.get_k());
        return word;
    }

    template <class Compare>
    static auto equal_range(const sdsl::int_vector<>& order, const Compare& compare) {
        auto bound = [&](bool upper) {
            size_t first = 0, last = order.size();
            while (first < last) {
                const auto middle = first + (last - first) / 2;
                const int comparison = compare(order[middle]);
                if (comparison < 0 || (upper && comparison == 0)) first = middle + 1;
                else last = middle;
            }
            return first;
        };
        return std::pair { bound(false), bound(true) };
    }

    template <class Dictionary>
    size_t call(const DBGSSHash& graph, const Dictionary& dictionary,
              std::string_view query, size_t minimum,
              const std::function<void(node_index, size_t)>& callback, size_t maximum) {
        using Kmer = typename Dictionary::kmer_type;
        const auto valid_end = std::find_if(query.begin(), query.end(),
                                           [](char c) { return !Kmer::is_valid(c); });
        query = query.substr(0, valid_end - query.begin());
        if (query.size() < minimum || !dictionary.num_kmers())
            return 0;
        std::call_once(initialized, [&] { build<true>(graph, dictionary, prefix_order); });

        Kmer wanted;
        for (size_t i = 0; i < query.size(); ++i)
            wanted.set(i, Kmer::char_to_uint(query[i]));
        auto compare = [&](node_index node) {
            const Kmer candidate = word_at(graph, dictionary, node);
            for (size_t i = 0; i < minimum; ++i) {
                if (candidate.at(i) < wanted.at(i)) return -1;
                if (candidate.at(i) > wanted.at(i)) return 1;
            }
            return 0;
        };
        const auto [first, last] = equal_range(prefix_order, compare);
        if (last - first > maximum)
            return last - first;
        for (size_t i = first; i < last; ++i) {
            const node_index node = prefix_order[i];
            const Kmer candidate = word_at(graph, dictionary, node);
            size_t length = minimum;
            while (length < query.size() && candidate.at(length) == wanted.at(length))
                ++length;
            callback(node, length);
        }
        return last - first;
    }

    template <class Dictionary>
    size_t call_suffix(const DBGSSHash& graph, const Dictionary& dictionary,
                     std::string_view query, size_t minimum,
                     const std::function<void(node_index, uint64_t)>& callback, size_t maximum) {
        using Kmer = typename Dictionary::kmer_type;
        if (!dictionary.num_kmers()
                || std::any_of(query.begin(), query.end(), [](char c) { return !Kmer::is_valid(c); }))
            return 0;
        std::call_once(suffix_initialized, [&] { build<false>(graph, dictionary, suffix_order); });

        Kmer wanted;
        for (size_t i = 0; i < query.size(); ++i)
            wanted.set(i, Kmer::char_to_uint(query[i]));
        for (size_t length = query.size(); length >= minimum; --length) {
            wanted.take_chars(length);
            auto compare = [&](node_index node) {
                auto candidate = word_at(graph, dictionary, node);
                candidate.drop_chars(graph.get_k() - length);
                return candidate.bits < wanted.bits ? -1 : candidate.bits > wanted.bits ? 1 : 0;
            };
            const auto [first, last] = equal_range(suffix_order, compare);
            if (first == last)
                continue;
            if (last - first > maximum)
                return maximum + 1;
            for (size_t i = first; i < last; ++i)
                callback(suffix_order[i], length);
            return last - first;
        }
        return 0;
    }
};

size_t DBGSSHash::call_nodes_with_prefix_matching(
        std::string_view query, size_t minimum,
        const std::function<void(node_index, size_t)>& callback, size_t maximum) const {
    query = query.substr(0, k_);
    minimum = std::max(size_t(1), minimum);
    if (query.size() < minimum)
        return 0;
    return std::visit([&](const auto& dictionary) {
        return seed_index_cache_->call(*this, dictionary, query, minimum, callback, maximum);
    }, dict_);
}

size_t DBGSSHash::call_nodes_with_suffix_matching_longest_prefix(
        std::string_view query, std::function<void(node_index, uint64_t)> callback,
        size_t minimum, size_t maximum) const {
    assert(query.size() <= k_);
    minimum = std::max(size_t(1), minimum);
    if (query.size() < minimum)
        return 0;
    return std::visit([&](const auto& dictionary) {
        return seed_index_cache_->call_suffix(*this, dictionary, query, minimum, callback, maximum);
    }, dict_);
}

size_t DBGSSHash::dict_size() const {
    return std::visit([](const auto& d) { return d.num_kmers(); }, dict_);
}

DBGSSHash::DBGSSHash(size_t k, Mode mode)
      : seed_index_cache_(std::make_shared<SeedIndexCache>()),
        k_(k), num_nodes_(0), mode_(mode) {
    if (mode == PRIMARY)
        mode_ = CANONICAL;

    size_t odd_k = (k_ | 1);

    if (odd_k * bits_per_char <= 64) {
        dict_.emplace<sshash::dictionary<kmer_t<KmerInt64>, sshash::decoded_offsets>>();
    } else if (odd_k * bits_per_char <= 128) {
        dict_.emplace<sshash::dictionary<kmer_t<KmerInt128>, sshash::decoded_offsets>>();
    } else if (odd_k * bits_per_char <= 256) {
        dict_.emplace<sshash::dictionary<kmer_t<KmerInt256>, sshash::decoded_offsets>>();
    } else {
        common::logger->error("k too big: {} > {}", odd_k, 256 / bits_per_char);
        throw std::length_error("k fail");
    }
}

std::pair<bit_vector_smart, bit_vector_smart> generate_succ_pred(const DBGSSHash& graph,
                                                                 uint64_t num_nodes) {
    sdsl::bit_vector succ_is_next(num_nodes + 1, false);
    sdsl::bit_vector pred_is_prev(num_nodes + 1, false);
    bool with_rc = (graph.get_mode() != DBGSSHash::BASIC);

    std::visit(
            [&](const auto& dict) {
                using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;

                ProgressBar progress_bar(num_nodes, "Checking succ/pred",
                                         std::cerr, !common::get_verbose());

                static const uint64_t BS = 1'048'576;

                // make sure each thread gets a disjoint block of words
                static_assert(BS % (sizeof(*succ_is_next.data()) * 8) == 0);

                std::atomic_thread_fence(std::memory_order_release);

#pragma omp parallel for num_threads(get_num_threads()) schedule(static)
                for (DBGSSHash::node_index begin = 0; begin <= num_nodes; begin += BS) {
                    DBGSSHash::node_index end = std::min(begin + BS, num_nodes + 1);

                    bool last_single_outdeg = false;
                    kmer_t kmer;
                    for (DBGSSHash::node_index node
                         = std::max(begin, DBGSSHash::node_index(1));
                         node < end; ++node) {
                        DBGSSHash::node_index next = DBGSSHash::npos;
                        if (!last_single_outdeg) {
                            std::string string_kmer = graph.get_node_sequence(node);
                            kmer = sshash::util::string_to_uint_kmer<kmer_t>(
                                    string_kmer.c_str(), graph.get_k());
                        }

                        last_single_outdeg = false;

                        auto nb = dict.kmer_neighbours(kmer, with_rc);
                        bool found_fw = false;
                        size_t found_i = 0;
                        for (size_t i = 0; i < nb.forward.size(); i++) {
                            if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                                if (found_fw) {
                                    next = DBGSSHash::npos;
                                    break;
                                } else {
                                    found_fw = true;
                                    next = DBGSSHash::sshash_to_graph_index(nb.forward[i].kmer_id);
                                    found_i = i;
                                    if (nb.forward[i].kmer_orientation < 0)
                                        next = graph.reverse_complement(next);
                                }
                            }
                        }

                        // The next index after the last stored node may denote
                        // a reverse-complement node, not the next stored k-mer.
                        if (node < num_nodes && next == node + 1) {
                            assert(found_fw);
                            succ_is_next[node] = true;
                            last_single_outdeg = true;
                            kmer.drop_char();
                            kmer.set(graph.get_k() - 1, found_i);
                        }

                        bool found_bw = false;
                        DBGSSHash::node_index prev = DBGSSHash::npos;
                        for (size_t i = 0; i < nb.backward.size(); i++) {
                            if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                                if (found_bw) {
                                    prev = DBGSSHash::npos;
                                    break;
                                } else {
                                    found_bw = true;
                                    prev = DBGSSHash::sshash_to_graph_index(nb.backward[i].kmer_id);
                                    if (nb.backward[i].kmer_orientation < 0)
                                        prev = graph.reverse_complement(prev);
                                }
                            }
                        }

                        if (prev != DBGSSHash::npos && prev == node - 1) {
                            assert(found_bw);
                            pred_is_prev[node] = true;
                        }
                    }

                    progress_bar += end - begin;
                }

                std::atomic_thread_fence(std::memory_order_acquire);
            },
            graph.data());

    return std::make_pair(bit_vector_smart(std::move(succ_is_next)),
                          bit_vector_smart(std::move(pred_is_prev)));
}

DBGSSHash::DBGSSHash(const std::string& input_filename,
                     size_t k,
                     Mode mode,
                     size_t num_chars,
                     const std::filesystem::path &swap_dir)
          : DBGSSHash(k, mode) {
    if (k <= 1)
        throw std::domain_error("k must be at least 2");

    if (mode == CANONICAL && (k % 2) == 0)
        throw std::domain_error("Primary graphs only supported for odd k");

    sshash::build_configuration build_config;
    build_config.k = k;

    // use a value of m recommended here:
    // https://twitter.com/giulio_pibiri/status/1803417277213114501
    // otherwise, use the next odd value greater than or equal to floor(k / 2)
    build_config.m
            = std::min({ (num_chars > 0
                                  ? uint64_t(ceil(log(static_cast<double>(num_chars))
                                                  / log(static_cast<double>(alphabet_.size())))
                                             + 1)
                                  : uint64_t(k / 2))
                                 | 1,
                         std::visit(
                                 [](const auto& dict) -> uint64_t {
                                     return std::decay_t<decltype(dict)>::kmer_type::max_m;
                                 },
                                 dict_),
                         uint64_t(k - 2) })
            | 1;

    build_config.verbose = common::get_verbose();
    build_config.num_threads = get_num_threads();
    // SSHash construction always uses temporary files. An empty caller path
    // selects the system temporary directory, just like the default argument.
    build_config.tmp_dirname = swap_dir.empty() ? std::filesystem::temp_directory_path()
                                               : swap_dir;
    build_config.allow_reverse_complement_duplicates = mode == BASIC;

    std::visit([&](auto& d) { d.build(input_filename, build_config); }, dict_);

    num_nodes_ = dict_size();

    std::tie(succ_is_next_, pred_is_prev_) = generate_succ_pred(*this, num_nodes_);
}

std::string DBGSSHash::file_extension() const {
    return kExtension;
}

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)>& on_insertion) {
    throw std::logic_error("adding sequences not supported");
}

sdsl::bit_vector DBGSSHash::get_invalid_char_mask(std::string_view sequence) const {
    sdsl::bit_vector invalid_char(sequence.size());

    std::visit([&](const auto& dict) {
        using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
        for (size_t i = 0; i < sequence.size(); ++i) {
            invalid_char[i] = !kmer_t::is_valid(sequence[i]);
        }
    }, dict_);

    return invalid_char;
}

template <bool with_rc, class Dict>
void map_to_nodes_with_rc_impl(const DBGSSHash& graph,
                               const Dict& dict,
                               std::string_view sequence,
                               const std::function<void(sshash::lookup_result)>& callback,
                               const std::function<bool()>& terminate) {
    size_t k = graph.get_k();
    size_t n = sequence.size();
    if (terminate() || n < k)
        return;

    if (!dict.num_kmers()) {
        for (size_t i = 0; i + k <= sequence.size() && !terminate(); ++i) {
            callback(sshash::lookup_result());
        }
        return;
    }

    using kmer_t = typename Dict::kmer_type;

    if (with_rc) {
        auto parser = sshash::streaming_query<Dict>(&dict);
        for (size_t i = 0; i + k <= sequence.size() && !terminate(); ++i) {
            auto ret_val = parser.lookup(sequence.data() + i);
            // Direct dictionary lookup assumes a valid alphabet. Streaming
            // lookup also accepts ambiguous bases and returns no match.
            assert(!sshash::util::is_valid<kmer_t>(sequence.data() + i, k)
                   || sshash::equal_lookup_result(ret_val,
                                               dict.lookup(sequence.data() + i, with_rc)));
            callback(ret_val);
        }
    } else {
        auto invalid_char = graph.get_invalid_char_mask(sequence);
        auto invalid_kmer = utils::drag_and_mark_segments(invalid_char, true, k);

        bool begin = true;
        kmer_t uint_kmer;
        for (size_t i = 0; i + k <= n && !terminate(); ++i) {
            sshash::lookup_result ret_val;
            if (invalid_kmer[i + k - 1]) {
                begin = true;
                callback(ret_val);
                continue;
            }

            if (begin) {
                begin = false;
                uint_kmer = sshash::util::string_to_uint_kmer<kmer_t>(sequence.data() + i, k);
            } else {
                uint_kmer.drop_char();
                uint_kmer.set(k - 1, kmer_t::char_to_uint(sequence[i + k - 1]));
            }
            ret_val = dict.lookup(uint_kmer, with_rc);

            assert(sshash::equal_lookup_result(ret_val,
                                               dict.lookup(sequence.data() + i, with_rc)));
            callback(ret_val);
        }
    }
}

template <bool with_rc>
void DBGSSHash::map_to_nodes_with_rc(std::string_view sequence,
                                     const std::function<void(node_index, bool)>& callback,
                                     const std::function<bool()>& terminate) const {
    std::visit(
            [&](const auto& dict) {
                map_to_nodes_with_rc_impl<with_rc>(
                        *this, dict, sequence,
                        [&](sshash::lookup_result res) {
                            auto node = sshash_to_graph_index(res.kmer_id);
                            assert(!node || in_graph(node));
                            callback(node, res.kmer_orientation < 0);
                        },
                        terminate);
            },
            dict_);
}

template void
DBGSSHash::map_to_nodes_with_rc<true>(std::string_view,
                                      const std::function<void(node_index, bool)>&,
                                      const std::function<bool()>&) const;
template void
DBGSSHash::map_to_nodes_with_rc<false>(std::string_view,
                                       const std::function<void(node_index, bool)>&,
                                       const std::function<bool()>&) const;

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)>& callback,
                             const std::function<bool()>& terminate) const {
    if (mode_ == BASIC) {
        map_to_nodes_sequentially(sequence, callback, terminate);
    } else {
        map_to_nodes_with_rc<true>(
                sequence, [&](node_index node, bool) { callback(node); }, terminate);
    }
}

DBGSSHash::node_index DBGSSHash::reverse_complement(node_index node) const {
    if (node == npos)
        return npos;

    assert(in_graph(node));
    if (node > dict_size())
        return node - dict_size();

    if (k_ % 2 == 1)
        return node + dict_size();

    std::string str_kmer(k_, ' ');
    uint64_t ssh_idx = graph_index_to_sshash(node);
    std::visit([&](const auto& d) { d.access(ssh_idx, str_kmer.data()); }, dict_);
    std::string rc_str_kmer(str_kmer);
    ::reverse_complement(rc_str_kmer.begin(), rc_str_kmer.end());
    return str_kmer != rc_str_kmer ? node + dict_size() : node;
}

DBGSSHash::node_index DBGSSHash::get_base_node(node_index node) const {
    if (node == npos)
        return npos;

    assert(in_graph(node));

    return node > dict_size() ? node - dict_size() : node;
}

uint64_t DBGSSHash::num_nodes() const {
    return mode_ != CANONICAL ? dict_size() : dict_size() * 2;
}

uint64_t DBGSSHash::max_index() const {
    return mode_ != BASIC ? dict_size() * 2 : dict_size();
}

void DBGSSHash::map_to_nodes_sequentially(std::string_view sequence,
                                          const std::function<void(node_index)>& callback,
                                          const std::function<bool()>& terminate) const {
    if (mode_ != BASIC) {
        map_to_nodes_with_rc<true>(
                sequence,
                [&](node_index n, bool revcompl) {
                    callback(n && revcompl ? reverse_complement(n) : n);
                },
                terminate);
    } else {
        map_to_nodes_with_rc<false>(
                sequence, [&](node_index node, bool) { callback(node); }, terminate);
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    assert(in_graph(node));
    // Fast path: if successor is the next index, just check if the character matches
    if (node < succ_is_next_.size()) {
        if (succ_is_next_[node]) {
            if (get_last_char(node + 1) == next_char)
                return node + 1;
            else
                return npos; // degree-1 node, but character doesn't match
        }
    } else {
        node_index rev_comp = reverse_complement(node);
        if (pred_is_prev_[rev_comp]) {
            --rev_comp;
            if (complement(get_first_char(rev_comp)) == next_char)
                return reverse_complement(rev_comp);
            else
                return npos;
        }
    }
    // Slow path: full MPHF lookup
    return kmer_to_node(get_node_sequence(node).substr(1) + next_char);
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));
    // TODO: if a node is in the middle of a unitig, then we only need to check the previous node index
    std::string string_kmer = prev_char + get_node_sequence(node);
    string_kmer.pop_back();
    return kmer_to_node(string_kmer);
}

template <bool with_rc>
void DBGSSHash::adjacent_outgoing_nodes_with_rc(
        node_index node,
        const std::function<void(node_index, bool)>& callback) const {
    assert(in_graph(node));

    if (node < succ_is_next_.size()) {
        if (succ_is_next_[node]) {
            callback(node + 1, false);
            return;
        }
    } else {
        node_index rev_comp = reverse_complement(node);
        if (pred_is_prev_[rev_comp--]) {
            callback(rev_comp, true);
            return;
        }
    }

    std::string kmer = get_node_sequence(node);
    std::visit(
            [&](const auto& dict) {
                auto nb = dict.kmer_forward_neighbours(kmer.c_str(), with_rc);
                for (size_t i = 0; i < nb.forward.size(); i++) {
                    if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                        callback(sshash_to_graph_index(nb.forward[i].kmer_id),
                                 nb.forward[i].kmer_orientation < 0);
                    }
                }
            },
            dict_);
}

template void DBGSSHash::adjacent_outgoing_nodes_with_rc<true>(
        node_index,
        const std::function<void(node_index, bool)>&) const;
template void DBGSSHash::adjacent_outgoing_nodes_with_rc<false>(
        node_index,
        const std::function<void(node_index, bool)>&) const;

template <bool with_rc>
void DBGSSHash::call_outgoing_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(in_graph(node));

    if (node < succ_is_next_.size()) {
        if (succ_is_next_[node]) {
            ++node;
            callback(node, get_last_char(node), false);
            return;
        }
    } else {
        node_index rev_comp = reverse_complement(node);
        if (pred_is_prev_[rev_comp--]) {
            callback(rev_comp, complement(get_first_char(rev_comp)), true);
            return;
        }
    }

    std::string kmer = get_node_sequence(node);
    std::visit(
            [&](const auto& dict) {
                using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
                auto nb = dict.kmer_forward_neighbours(kmer.c_str(), with_rc);
                for (size_t i = 0; i < nb.forward.size(); i++) {
                    if (nb.forward[i].kmer_id != sshash::constants::invalid_uint64) {
                        callback(sshash_to_graph_index(nb.forward[i].kmer_id),
                                 kmer_t::uint64_to_char(i), nb.forward[i].kmer_orientation < 0);
                    }
                }
            },
            dict_);
}

template void DBGSSHash::call_outgoing_kmers_with_rc<true>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;
template void DBGSSHash::call_outgoing_kmers_with_rc<false>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;


template <bool with_rc>
void DBGSSHash::adjacent_incoming_nodes_with_rc(
        node_index node,
        const std::function<void(node_index, bool)>& callback) const {
    assert(in_graph(node));

    if (node < pred_is_prev_.size()) {
        if (pred_is_prev_[node]) {
            callback(node - 1, false);
            return;
        }
    } else {
        node_index rev_comp = reverse_complement(node);
        if (succ_is_next_[rev_comp++]) {
            callback(rev_comp, true);
            return;
        }
    }

    std::string kmer = get_node_sequence(node);
    std::visit(
            [&](const auto& dict) {
                auto nb = dict.kmer_backward_neighbours(kmer.c_str(), with_rc);
                for (size_t i = 0; i < nb.backward.size(); i++) {
                    if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                        callback(sshash_to_graph_index(nb.backward[i].kmer_id),
                                 nb.backward[i].kmer_orientation < 0);
                    }
                }
            },
            dict_);
}

template void DBGSSHash::adjacent_incoming_nodes_with_rc<true>(
        node_index,
        const std::function<void(node_index, bool)>&) const;
template void DBGSSHash::adjacent_incoming_nodes_with_rc<false>(
        node_index,
        const std::function<void(node_index, bool)>&) const;


template <bool with_rc>
void DBGSSHash::call_incoming_kmers_with_rc(
        node_index node,
        const std::function<void(node_index, char, bool)>& callback) const {
    assert(in_graph(node));

    if (node < pred_is_prev_.size()) {
        if (pred_is_prev_[node]) {
            --node;
            callback(node, get_first_char(node), false);
            return;
        }
    } else {
        node_index rev_comp = reverse_complement(node);
        if (succ_is_next_[rev_comp++]) {
            callback(rev_comp, complement(get_last_char(rev_comp)), true);
            return;
        }
    }

    std::string kmer = get_node_sequence(node);
    std::visit(
            [&](const auto& dict) {
                using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
                auto nb = dict.kmer_backward_neighbours(kmer.c_str(), with_rc);
                for (size_t i = 0; i < nb.backward.size(); i++) {
                    if (nb.backward[i].kmer_id != sshash::constants::invalid_uint64) {
                        callback(sshash_to_graph_index(nb.backward[i].kmer_id),
                                 kmer_t::uint64_to_char(i), nb.backward[i].kmer_orientation < 0);
                    }
                }
            },
            dict_);
}


template void DBGSSHash::call_incoming_kmers_with_rc<true>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;
template void DBGSSHash::call_incoming_kmers_with_rc<false>(
        node_index,
        const std::function<void(node_index, char, bool)>&) const;

void DBGSSHash::call_outgoing_kmers(node_index node,
                                    const OutgoingEdgeCallback& callback) const {
    if (mode_ != BASIC) {
        call_outgoing_kmers_with_rc<true>(node, [&](node_index next, char c, bool revcompl) {
            callback(revcompl ? reverse_complement(next) : next, c);
        });
    } else {
        call_outgoing_kmers_with_rc<false>(node, [&](node_index next, char c, bool) {
            callback(next, c);
        });
    }
}

void DBGSSHash::adjacent_outgoing_nodes(node_index node,
                                        const std::function<void(node_index)>& callback) const {
    if (mode_ != BASIC) {
        adjacent_outgoing_nodes_with_rc<true>(node, [&](node_index next, bool revcompl) {
            callback(revcompl ? reverse_complement(next) : next);
        });
    } else {
        adjacent_outgoing_nodes_with_rc<false>(node, [&](node_index next, bool) {
            callback(next);
        });
    }
}

void DBGSSHash::call_incoming_kmers(node_index node,
                                    const IncomingEdgeCallback& callback) const {
    if (mode_ != BASIC) {
        call_incoming_kmers_with_rc<true>(node, [&](node_index prev, char c, bool revcompl) {
            callback(revcompl ? reverse_complement(prev) : prev, c);
        });
    } else {
        call_incoming_kmers_with_rc<false>(node, [&](node_index prev, char c, bool) {
            callback(prev, c);
        });
    }
}

void DBGSSHash::adjacent_incoming_nodes(node_index node,
                                        const std::function<void(node_index)>& callback) const {
    if (mode_ != BASIC) {
        adjacent_incoming_nodes_with_rc<true>(node, [&](node_index prev, bool revcompl) {
            callback(revcompl ? reverse_complement(prev) : prev);
        });
    } else {
        adjacent_incoming_nodes_with_rc<false>(node, [&](node_index prev, bool) {
            callback(prev);
        });
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    assert(in_graph(node));
    size_t res = 0;
    adjacent_outgoing_nodes(node, [&](node_index) { ++res; });
    return res;
}

size_t DBGSSHash::indegree(node_index node) const {
    assert(in_graph(node));
    size_t res = 0;
    adjacent_incoming_nodes(node, [&](node_index) { ++res; });
    return res;
}

void DBGSSHash::call_nodes(const std::function<void(node_index)>& callback,
                           const std::function<bool()>& terminate,
                           size_t num_threads,
                           size_t batch_size) const {
#pragma omp parallel for num_threads(num_threads) schedule(static, batch_size)
    for (size_t node_idx = 1; node_idx <= dict_size(); ++node_idx) {
        if (terminate())
            continue;

        callback(node_idx);
    }

    if (terminate())
        return;

    if (mode_ == CANONICAL) {
#pragma omp parallel for num_threads(num_threads) schedule(static, batch_size)
        for (size_t node_idx = 1; node_idx <= dict_size(); ++node_idx) {
            if (terminate())
                continue;

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

    return std::visit(
            [&](const auto& d) {
                auto res = d.lookup(kmer.data(), with_rc);
                return std::make_pair(sshash_to_graph_index(res.kmer_id),
                                      res.kmer_orientation < 0);
            },
            dict_);
}
template std::pair<DBGSSHash::node_index, bool>
        DBGSSHash::kmer_to_node_with_rc<true>(std::string_view) const;
template std::pair<DBGSSHash::node_index, bool>
        DBGSSHash::kmer_to_node_with_rc<false>(std::string_view) const;

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    if (mode_ == CANONICAL) {
        auto res = kmer_to_node_with_rc<true>(kmer);
        return res.second ? reverse_complement(res.first) : res.first;
    } else {
        return kmer_to_node_with_rc<false>(kmer).first;
    }
}

char DBGSSHash::get_last_char(node_index node) const {
    assert(in_graph(node));
    assert(node <= dict_size());
    uint64_t ssh_idx = graph_index_to_sshash(node);
    return std::visit(
            [&](const auto& d) {
                using kmer_t = typename std::decay_t<decltype(d)>::kmer_type;
                auto [_, offset] = d.strings_offsets().id_to_offset(ssh_idx, get_k());
                sshash::kmer_iterator<kmer_t, bits::bit_vector> kmer_it(
                        d.strings(), get_k(),
                        kmer_t::bits_per_char * (offset + get_k() - 1));
                return kmer_t::uint64_to_char(kmer_it.get_next_char());
            },
            dict_);
}

char DBGSSHash::get_first_char(node_index node) const {
    assert(in_graph(node));
    assert(node <= dict_size());
    uint64_t ssh_idx = graph_index_to_sshash(node);
    return std::visit(
            [&](const auto& d) {
                using kmer_t = typename std::decay_t<decltype(d)>::kmer_type;
                auto [_, offset] = d.strings_offsets().id_to_offset(ssh_idx, get_k());
                sshash::kmer_iterator<kmer_t, bits::bit_vector> kmer_it(d.strings(), get_k(),
                                                       kmer_t::bits_per_char * offset);
                return kmer_t::uint64_to_char(kmer_it.get_next_char());
            },
            dict_);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    assert(in_graph(node));
    std::string str_kmer(k_, ' ');
    node_index node_canonical = node > dict_size() ? node - dict_size() : node;
    uint64_t ssh_idx = graph_index_to_sshash(node_canonical);
    std::visit([&](const auto& d) { d.access(ssh_idx, str_kmer.data()); }, dict_);

    if (node > node_canonical)
        ::reverse_complement(str_kmer.begin(), str_kmer.end());

    return str_kmer;
}

std::vector<DBGSSHash::seq_index> DBGSSHash::get_sequence_ids(node_index node) const {
    assert(node <= max_index());
    if (node == npos)
        return std::vector<seq_index>(1, node);

    if (node > dict_size()) {
        assert(mode_ == CANONICAL);
        size_t num_base_contigs = succ_is_next_.size() + pred_is_prev_.size() - succ_is_next_.num_set_bits() - pred_is_prev_.num_set_bits();
        auto ret_val = get_sequence_ids(reverse_complement(node));
        for (seq_index &id : ret_val) {
            id += num_base_contigs;
        }

        return ret_val;
    }

    // treat each unitig as a sequence
    size_t srank = succ_is_next_.rank0(node);
    size_t prank = pred_is_prev_.rank0(node);

    seq_index contig_id = srank + prank;

    bool s_unitig_back = !succ_is_next_[node];
    bool s_unitig_front = !pred_is_prev_[node];

    if (s_unitig_back && !s_unitig_front) {
        assert(contig_id);
        --contig_id;
    }

    return std::vector<seq_index>(1, contig_id + 1);
}

void DBGSSHash::serialize(std::ostream& out) const {
    essentials::generic_saver saver(out);

    saver.visit(num_nodes_);
    saver.visit(k_);
    saver.visit(mode_);

    if (num_nodes()) {
        std::visit([&](const auto& d) { saver.visit(d); }, dict_);
        succ_is_next_.serialize(out);
        pred_is_prev_.serialize(out);
    }
}

void DBGSSHash::serialize(const std::string& filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    std::ofstream fout(suffixed_filename, std::ios::binary);

    if (!fout.good()) {
        common::logger->error("Failed to create file {}", suffixed_filename);
        throw std::runtime_error("Invalid out file");
    }

    serialize(fout);
}

bool DBGSSHash::load(std::istream& in) {
    if (!in.good())
        return false;

    struct RestoreExceptions {
        std::istream &stream;
        std::ios::iostate mask;
        ~RestoreExceptions() {
            try {
                stream.exceptions(mask);
            } catch (const std::ios_base::failure &) {
                // The mask is restored even if the failed read still sets a flag.
            }
        }
    } restore_exceptions { in, in.exceptions() };

    try {
        // Stop on a short read before the loader uses an unread array length.
        in.exceptions(in.badbit | in.failbit | in.eofbit);
        essentials::generic_loader loader(in);
        size_t num_nodes;
        size_t k;
        Mode mode;
        loader.visit(num_nodes);
        loader.visit(k);
        loader.visit(mode);

        if (mode == PRIMARY) {
            common::logger->warn("Graph constructed with old version and is incompatible. Please rebuild only graph (not annotations).");
            throw std::runtime_error("Incompatible graph");
        }

        DBGSSHash loaded(k, mode);
        loaded.num_nodes_ = num_nodes;

        if (num_nodes) {
            // Enable zero-copy mmap if the stream is backed by an mmap'd file
            if (auto *mmap_stream = dynamic_cast<sdsl::mmap_ifstream*>(&in)) {
                auto ctx = mmap_stream->get_mmap_context();
                loader.set_mmap(ctx->data(), ctx->file_size_bytes(), ctx);
            }
            std::visit([&](auto& d) { d.visit(loader); }, loaded.dict_);
            if (!loaded.succ_is_next_.load(in) || !loaded.pred_is_prev_.load(in)
                    || loaded.dict_size() != num_nodes
                    || loaded.succ_is_next_.size() != num_nodes + 1
                    || loaded.pred_is_prev_.size() != num_nodes + 1)
                return false;

            // Earlier v6 graphs omitted the final predecessor cache bit.
            // Repair that bit without repeating the graph-wide neighbor scan.
            if (num_nodes > 1 && !loaded.pred_is_prev_[num_nodes]) {
                node_index predecessor = npos;
                size_t num_predecessors = 0;
                loaded.adjacent_incoming_nodes(num_nodes, [&](node_index node) {
                    predecessor = node;
                    ++num_predecessors;
                });
                if (num_predecessors == 1 && predecessor == num_nodes - 1) {
                    auto predecessors = loaded.pred_is_prev_.to_vector();
                    predecessors[num_nodes] = true;
                    loaded.pred_is_prev_ = bit_vector_smart(std::move(predecessors));
                }
            }
        }

        if (!in.good())
            return false;

        *this = std::move(loaded);
        return true;
    } catch (const std::exception &error) {
        common::logger->error("Failed to load SSHash graph: {}", error.what());
        return false;
    }
}

bool DBGSSHash::load(const std::string& filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    auto fin = utils::open_ifstream(suffixed_filename);
    return load(*fin);
}

} // namespace graph
} // namespace mtg
