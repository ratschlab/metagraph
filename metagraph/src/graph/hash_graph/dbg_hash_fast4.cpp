#include "dbg_hash_fast4.hpp"

#include <cassert>

#include <tsl/ordered_set.h>
//#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <libmaus2/util/NumberSerialisation.hpp>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"

#include <unordered_set> //TODO remove after debug

template <typename KMER = KmerExtractor2Bit::Kmer64>
class DBGHashFast4Impl : public DBGHashFast4::DBGHashFast4Interface {
    using Kmer = KMER;
    using TAlphabet = KmerExtractor2Bit::TAlphabet;
    using KmerPrefix = typename KMER::WordType;

#if _PROTEIN_GRAPH
    using Bits = uint32_t;
#elif _DNA_CASE_SENSITIVE_GRAPH
    using Bits = uint16_t;
#elif _DNA5_GRAPH
    using Bits = unsigned char;
#elif _DNA_GRAPH
    using Bits = unsigned char;
#else
    static_assert(false, "invalid or undefined alphabet");
#endif

    //using KmerIndex = tsl::robin_map<KmerPrefix,
    //                                 Bits,
    //                                 utils::Hash<KmerPrefix>,
    //                                 std::equal_to<KmerPrefix>,
    //                                 std::allocator<std::pair<KmerPrefix, Bits>>,
    //                                 false,
    //                                 tsl::rh::power_of_two_growth_policy<2>>;
    using KmerIndex = tsl::ordered_set<KmerPrefix,
                                     utils::Hash<KmerPrefix>,
                                     std::equal_to<KmerPrefix>,
                                     std::allocator<KmerPrefix>,
                                     //std::deque<KmerPrefix, std::allocator<KmerPrefix>>,
                                     std::vector<KmerPrefix, std::allocator<KmerPrefix>>,
                                     std::uint64_t>;

    using KmerIterator = typename KmerIndex::iterator;
    using KmerConstIterator = typename KmerIndex::const_iterator;

  public:
    explicit DBGHashFast4Impl(size_t k,
                              bool canonical_mode,
                              bool packed_serialization);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const;
    void call_nodes(const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &stop_early) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index) const;
    bool has_multiple_outgoing(node_index) const;

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index) const;

    node_index kmer_to_node(const std::string &kmer) const;

    std::string get_node_sequence(node_index node) const;

    size_t get_k() const { return k_; }
    bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const { return kmers_.size() << (KMER::kBitsPerChar * 4); }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const;

    const std::string& alphabet() const { return seq_encoder_.alphabet; }

    bool in_graph(node_index node) const;

  private:
    Vector<Kmer> sequence_to_kmers(const std::string &sequence, bool canonical = false) const {
        return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, canonical);
    }

    bool kmers_overlap(const Kmer &out_kmer, const Kmer &in_kmer) const {
        KMER overlap = out_kmer;
        overlap.to_next(k_, in_kmer[k_ - 1]);
        return overlap == in_kmer;
    }

    size_t find_shift(node_index node) const {
        return ((node - 1) >> (KMER::kBitsPerChar * 3)) & ((1 << KMER::kBitsPerChar) - 1);
    }

    size_t get_rest(KMER kmer, size_t shift) const {
        const KmerPrefix before_mask = before_mask_[shift];
        const KmerPrefix after_mask = after_mask_[shift];

        return ((KmerPrefix(kmer.data()) & before_mask) << int(KMER::kBitsPerChar * (3 - shift))) |
            ((KmerPrefix(kmer.data()) & after_mask) >> int(KMER::kBitsPerChar * (k_ - 3 + shift)));
    }

    void print_internal_representation() {
        for (auto iter = kmers_.begin(); iter != kmers_.end(); ++iter) {
            auto val = values_[iter - kmers_.begin()];
            std::cout << KMER(*iter).to_string(k_ - 3, seq_encoder_.alphabet) << " " << std::bitset<64>(val[0]) << "," << std::bitset<64>(val[1]) << "," << std::bitset<64>(val[2]) << "," << std::bitset<64>(val[3]) << " ";
            for (size_t shift = 0; shift <= 3; ++shift) {
                for (size_t i = 0; i < (1ul << (KMER::kBitsPerChar * 3)); ++i) {
                    if (val[shift] & (1ull << i))
                        std::cout << KMER(i).to_string(3, seq_encoder_.alphabet) << " (" << shift << "),";
                }
            }
            std::cout << std::endl;
        }
        std::cout << "---" << std::endl;
    }

    node_index get_index(const KmerConstIterator &iter) const;
    node_index get_index(const Kmer &kmer) const;
    node_index get_index(const Kmer &kmer, const KmerConstIterator &iter) const;

    KmerIterator get_iter(node_index node);
    KmerConstIterator get_const_iter(node_index node) const;

    const Kmer get_kmer(node_index node,
                        KmerConstIterator kmer_iter) const;

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;
    std::vector<Bits> bits_;
    KmerExtractor2Bit seq_encoder_;

    bool packed_serialization_;

    const KmerPrefix kIgnoreLastCharMask;
    const Bits kLastCharMask = (Bits(1) << seq_encoder_.alphabet.size()) - 1;
    const Bits kIncomingEdgesMask = kLastCharMask << seq_encoder_.alphabet.size();
    static constexpr auto kExtension = DBGHashFast4::kExtension;

    unsigned char prefixes_[4096];
    std::vector<std::array<uint64_t, 4>> values_;
    std::vector<bool> may_contain_source_kmer_;

    const KmerPrefix prefix_mask_ = KmerPrefix((1 << KMER::kBitsPerChar * 6) - 1);
    const KmerPrefix first3_mask_ = prefix_mask_ >> int(KMER::kBitsPerChar * 3);
    KmerPrefix before_mask_[4];
    KmerPrefix after_mask_[4];
    KmerPrefix both_mask_[4];
    KmerPrefix neither_mask_[4];

    std::unordered_set<std::string> debug_;
};

template <typename KMER>
DBGHashFast4Impl<KMER>::DBGHashFast4Impl(size_t k,
                                         bool canonical_mode,
                                         bool packed_serialization)
      : k_(k),
        canonical_mode_(canonical_mode),
        packed_serialization_(packed_serialization),
        kIgnoreLastCharMask(~(KmerPrefix((1llu << KMER::kBitsPerChar) - 1) << int(KMER::kBitsPerChar * (k - 1))))
{
    for (uint32_t i = 0; i < (KmerPrefix(1) << (KMER::kBitsPerChar * 6)); ++i) {
        KMER prefix{KmerPrefix(i)};

        const KmerPrefix lower_2 = (1 << 3) - 1;
        
        for (unsigned char shift = 0; shift <= 4; ++shift) {
            const int shift_bits = shift * KMER::kBitsPerChar;
            KmerPrefix data = prefix.data();

            if (shift == 4 ||
                (__builtin_popcount(data & (lower_2 << shift_bits)) + 
                 __builtin_popcount(data & (lower_2 << (shift_bits + 3)))) % 3 == 1) {
                if (shift == 4) {
                    if ((__builtin_popcount(data & (lower_2 << (1 * KMER::kBitsPerChar))) + 
                         __builtin_popcount(data & (lower_2 << (2 * KMER::kBitsPerChar)))) % 2 == 1) {
                        prefixes_[size_t(data)] = 2;
                    } else if ((__builtin_popcount(data & (lower_2 << (0 * KMER::kBitsPerChar))) + 
                                __builtin_popcount(data & (lower_2 << (1 * KMER::kBitsPerChar)))) % 2 == 1) {
                        prefixes_[size_t(data)] = 1;
                    } else {
                        prefixes_[size_t(data)] = 0;
                    }
                } else {
                    prefixes_[size_t(data)] = shift;
                }
                break;
            }
        }
    }

    for (size_t shift = 0; shift <= 3; ++shift) {
        before_mask_[shift] = first3_mask_ >> int(KMER::kBitsPerChar * (3 - shift));
        after_mask_[shift] = (first3_mask_ >> int(KMER::kBitsPerChar * shift)) <<
                             int(KMER::kBitsPerChar * (k_ - (3 - shift)));
        both_mask_[shift] = before_mask_[shift] | after_mask_[shift];
        neither_mask_[shift] = ~both_mask_[shift];
    }
};


template <typename KMER>
void DBGHashFast4Impl<KMER>::add_sequence(const std::string &sequence,
                                            bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);
    kmers_.reserve(50'000'000);
    values_.reserve(50'000'000);
    may_contain_source_kmer_.reserve(50'000'000);

    size_t shift = 4;
    size_t last_shift = 0;
    KmerIterator iter = kmers_.end();
    auto kmers = sequence_to_kmers(sequence);
    for (size_t kmer_index = 0; kmer_index < kmers.size(); ++kmer_index) {
        const auto &kmer = kmers[kmer_index];
        bool inserted = false;

        shift = prefixes_[size_t(KmerPrefix(kmer.data()) & prefix_mask_)];
        if (shift + 1 != last_shift) {
            iter = kmers_.end();
        } else {
        }
        last_shift = shift;

        const size_t rest = get_rest(kmer, shift);

        bool may_contain_source_kmer = false;
        if (kmer_index > 0 && !kmers_overlap(kmers[kmer_index - 1], kmer)) {
            iter = kmers_.end();
            may_contain_source_kmer = true;
        }
        if (kmer_index == 0)
            may_contain_source_kmer = true;

        if (iter == kmers_.end()) {
            KmerPrefix key = KmerPrefix(kmer.data()) & neither_mask_[shift];
            key = key >> int(KMER::kBitsPerChar * shift);
            std::tie(iter, inserted) = kmers_.insert(key);
        }

        const uint64_t val = 1ull << rest;
        if (inserted) {
            std::array<uint64_t, 4> arr = {0, 0, 0, 0};
            arr[shift] = val;
            values_.push_back(arr);
            may_contain_source_kmer_.push_back(may_contain_source_kmer);
        } else {
            values_[iter - kmers_.begin()][shift] |= val;
            if (may_contain_source_kmer)
                may_contain_source_kmer_[iter - kmers_.begin()] = may_contain_source_kmer;
        }

/*
        print_internal_representation();
        call_nodes([&](node_index node) {
            //auto kmer_iter = get_iter(node);
            //KmerPrefix key = kmer_iter.key();
            //KmerPrefix rest = KmerPrefix(node - 1) & first3_mask_;
            std::cout << node << ":" << get_node_sequence(node) << " (" << find_shift(node) << ")" << ", ";
        }, [](){ return false; });
        std::cout << std::endl;
*/

        assert(kmer.to_string(k_, seq_encoder_.alphabet) == get_node_sequence(get_index(iter) + (shift << (KMER::kBitsPerChar * 3)) + rest));
        assert((get_index(iter) + (shift << (KMER::kBitsPerChar * 3)) + rest) == kmer_to_node(get_node_sequence(get_index(iter) + (shift << (KMER::kBitsPerChar * 3)) + rest)));
    }

    std::ignore = nodes_inserted;

/*
    auto kmers = sequence_to_kmers(sequence);
    KmerIterator next_kmer_prefix_it = kmers_.end();
    for (auto kmer_it = kmers.rbegin(); kmer_it != kmers.rend(); ++kmer_it) {
        const auto &kmer = *kmer_it;

        if (kmer_it == kmers.rbegin() || !kmers_overlap(kmer, *(kmer_it - 1))) {
            const auto &next_kmer_prefix = KmerPrefix(kmer.data()) >> KMER::kBitsPerChar;
            next_kmer_prefix_it = find_kmer(next_kmer_prefix);

            if (next_kmer_prefix_it == kmers_.end()) {
                next_kmer_prefix_it = kmers_.insert(next_kmer_prefix).first; // dummy
                bits_.push_back(Bits(0));
            }
        }

        Bits incoming_val = (Bits(1) << kmer[0]) << seq_encoder_.alphabet.size();
        assert(next_kmer_prefix_it != kmers_.end());
        bits_[get_bit_index(next_kmer_prefix_it)] |= incoming_val;

        Bits val = Bits(1) << kmer[k_ - 1];

        const auto &&key = KmerPrefix(kmer.data()) & kIgnoreLastCharMask;

        auto &&[iter, inserted] = kmers_.insert(key);
        assert(iter != kmers_.end());
        assert(iter == find_kmer(key));

        if (inserted)
            bits_.push_back(val);
        else
            bits_[get_bit_index(iter)] |= val;

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);

        next_kmer_prefix_it = iter;
    }

    if (!canonical_mode_)
        return;

    auto rev_kmers = sequence_to_kmers(seq_encoder_.reverse_complement(sequence));
    next_kmer_prefix_it = kmers_.end();
    for (auto kmer_it = rev_kmers.rbegin(); kmer_it != rev_kmers.rend(); ++kmer_it) {
        const auto &kmer = *kmer_it;

        if (kmer_it == kmers.rbegin() || !kmers_overlap(kmer, *(kmer_it - 1))) {
            const auto &next_kmer_prefix = KmerPrefix(kmer.data()) >> KMER::kBitsPerChar;
            next_kmer_prefix_it = find_kmer(next_kmer_prefix);

            if (next_kmer_prefix_it == kmers_.end()) {
                next_kmer_prefix_it = kmers_.insert(next_kmer_prefix).first; // dummy
                bits_.push_back(Bits(0));
            }
        }

        Bits incoming_val = (Bits(1) << kmer[0]) << seq_encoder_.alphabet.size();
        assert(next_kmer_prefix_it != kmers_.end());
        bits_[get_bit_index(next_kmer_prefix_it)] |= incoming_val;

        Bits val = Bits(1) << kmer[k_ - 1];

        const auto &&key = KmerPrefix(kmer.data()) & kIgnoreLastCharMask;

        auto &&[iter, inserted] = kmers_.insert(key);
        assert(iter != kmers_.end());
        assert(iter == find_kmer(key));

        if (inserted)
            bits_.push_back(val);
        else
            bits_[get_bit_index(iter)] |= val;

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);

        next_kmer_prefix_it = iter;
    }
*/
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
template <typename KMER>
void DBGHashFast4Impl<KMER>::map_to_nodes_sequentially(
                              std::string::const_iterator begin,
                              std::string::const_iterator end,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate) const {
    std::string sequence(begin, end);

    const auto &kmers = sequence_to_kmers(sequence);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {

        assert(it != kmers.end() || !is_valid);

        if (terminate())
            return;

        //assert((find_kmer(KmerPrefix((*it).data()) & kIgnoreLastCharMask) == kmers_.end())
        //       || get_index(*it) == npos
        //       || *it == get_kmer(get_index(*it), get_const_iter(get_index(*it))));

        callback(is_valid ? get_index(*it++) : npos);
    }
    assert(it == kmers.end());
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
template <typename KMER>
void DBGHashFast4Impl<KMER>::map_to_nodes(const std::string &sequence,
                                            const std::function<void(node_index)> &callback,
                                            const std::function<bool()> &terminate) const {

    const auto &kmers = sequence_to_kmers(sequence, canonical_mode_);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {

        assert(it != kmers.end() || !is_valid);

        if (terminate())
            return;

        assert(!get_index(*it) || *it == get_kmer(get_index(*it), get_const_iter(get_index(*it))));

        callback(is_valid ? get_index(*it++) : npos);
    }
    assert(it == kmers.end());
}

template <typename KMER>
void DBGHashFast4Impl<KMER>::call_outgoing_kmers(node_index node,
                                                 const OutgoingEdgeCallback &callback) const {
    KmerConstIterator iter = get_const_iter(node);
    KmerPrefix key = iter.key();
    std::ignore = key;
    KmerPrefix full = get_kmer(node, get_const_iter(node)).data();

    Kmer next_kmer = KMER(full);
    next_kmer.to_next(k_, 0);

    size_t next_shift = prefixes_[size_t(KmerPrefix(next_kmer.data()) & prefix_mask_)];

    KmerPrefix next_key = KmerPrefix(next_kmer.data()) & neither_mask_[next_shift];
    next_key = next_key >> int(KMER::kBitsPerChar * next_shift);

    assert(prefixes_[size_t(full & prefix_mask_)] == find_shift(node));
    if (next_shift + 1 == find_shift(node)) {
        assert(next_key == key);

        uint64_t val = values_[iter - kmers_.begin()][next_shift];

        for (size_t i = 0; i < seq_encoder_.alphabet.size(); ++i) {
            Kmer next_kmer = KMER(full);
            next_kmer.to_next(k_, i);

            const size_t next_rest = get_rest(next_kmer, next_shift);

            if (val & (1ull << next_rest)) {
                node_index base_node_index = get_index(iter);
                callback(base_node_index + (next_shift << (KMER::kBitsPerChar * 3)) + next_rest, seq_encoder_.decode(i));
            }
        }
    } else {
        KmerConstIterator next_iter;
        if (next_key == key) {
            next_iter = iter;
        } else {
            next_iter = kmers_.find(next_key);
        }

        for (size_t i = 0; i < seq_encoder_.alphabet.size(); ++i) {
            Kmer next_kmer = KMER(full);
            next_kmer.to_next(k_, i);
            next_shift = prefixes_[size_t(KmerPrefix(next_kmer.data()) & prefix_mask_)];

            KmerPrefix next_key_i = KmerPrefix(next_kmer.data()) & neither_mask_[next_shift];
            next_key_i = next_key_i >> int(KMER::kBitsPerChar * next_shift);

            node_index next_node;
            if (next_key_i == next_key) {
                next_node = get_index(next_kmer, next_iter);
                assert(get_index(next_kmer) == next_node);
            } else {
                next_node = get_index(next_kmer);
            }

            if (next_node != npos)
                callback(next_node, seq_encoder_.decode(i));
        }
    }
}

template <typename KMER>
void DBGHashFast4Impl<KMER>::call_incoming_kmers(node_index node,
                                                 const IncomingEdgeCallback &callback) const {
    KmerConstIterator iter = get_const_iter(node);
    KmerPrefix key = iter.key();
    KmerPrefix full = get_kmer(node, get_const_iter(node)).data();

    for (size_t i = 0; i < seq_encoder_.alphabet.size(); ++i) {
        Kmer prev_kmer = KMER(full);
        prev_kmer.to_prev(k_, i);
            
        size_t prev_shift = prefixes_[size_t(KmerPrefix(prev_kmer.data()) & prefix_mask_)];

        KmerPrefix prev_key = KmerPrefix(prev_kmer.data()) & neither_mask_[prev_shift];
        prev_key = prev_key >> int(KMER::kBitsPerChar * prev_shift);

        if (prev_key == key) {
            const size_t prev_rest = get_rest(prev_kmer, prev_shift);

            uint64_t val = values_[iter - kmers_.begin()][prev_shift];

            if (val & (1ull << prev_rest)) {
                node_index base_node_index = get_index(iter);
                callback(base_node_index + (prev_shift << (KMER::kBitsPerChar * 3)) + prev_rest, seq_encoder_.decode(i));
            }
        } else {
            node_index prev_node = get_index(prev_kmer);
            if (prev_node != npos)
                callback(prev_node, seq_encoder_.decode(i));
        }
    }
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::traverse(node_index node, char next_char) const {
    assert(in_graph(node));
    node_index ret = npos;
    call_outgoing_kmers(
        node,
        [&](node_index out_node, char c) {
        if (c == next_char)
            ret = out_node;
        }
    );
    return ret;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));
    node_index ret = npos;
    call_incoming_kmers(
        node,
        [&](node_index in_node, char c) {
        if (c == prev_char)
            ret = in_node;
        }
    );
    return ret;
}

template <typename KMER>
void
DBGHashFast4Impl<KMER>
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

template <typename KMER>
void
DBGHashFast4Impl<KMER>
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

template <typename KMER>
size_t DBGHashFast4Impl<KMER>::outdegree(node_index node) const {
    assert(in_graph(node));

    size_t outdegree = 0;
    call_outgoing_kmers(node, [&](auto, auto) { outdegree++; });
    return outdegree;
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) == 1;
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) > 1;
}

template <typename KMER>
size_t DBGHashFast4Impl<KMER>::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;
    call_incoming_kmers(node, [&](auto, auto) { indegree++; });
    return indegree;
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    if (!may_contain_source_kmer_[get_const_iter(node) - kmers_.begin()])
        return false;

    return indegree(node) == 0;
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    return indegree(node) == 1;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    return get_index(seq_encoder_.encode(kmer));
}

template <typename KMER>
std::string DBGHashFast4Impl<KMER>::get_node_sequence(node_index node) const {
    assert(in_graph(node));

    return seq_encoder_.kmer_to_sequence(get_kmer(node, get_const_iter(node)), k_);
}

class Serializer {
  public:
    explicit Serializer(std::ostream &os) : os_(os) {}

    template <class T>
    void operator()(const T &value) {
        os_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

  private:
    std::ostream &os_;
};

class Deserializer {
  public:
    explicit Deserializer(std::istream &is) : is_(is) {}

    template <class T>
    T operator()() {
        T value;
        is_.read(reinterpret_cast<char *>(&value), sizeof(T));
        return value;
    }

  private:
    std::istream &is_;
};

template <typename KMER>
void DBGHashFast4Impl<KMER>::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    out.exceptions(out.badbit | out.failbit);

    serialize_number(out, k_);

    Serializer serializer(out);

    if (packed_serialization_) {
        serialize_number(out, kmers_.size());
        std::for_each(kmers_.begin(), kmers_.end(), serializer);
        serialize_number(out, values_.size());
        std::for_each(values_.begin(), values_.end(), serializer);
        serialize_number(out, may_contain_source_kmer_.size());
        std::for_each(may_contain_source_kmer_.begin(), may_contain_source_kmer_.end(), serializer);
        //serialize_number(out, shifts_.size());
        //std::for_each(shifts_.begin(), shifts_.end(), serializer);
    } else {
        throw std::runtime_error("not implemented");
        //serialize_number(out, std::numeric_limits<uint64_t>::max());
        //kmers_.serialize(serializer);
    }

    serialize_number(out, canonical_mode_);
}

template <typename KMER>
void DBGHashFast4Impl<KMER>::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::load(std::istream &in) {
    if (!in.good())
        return false;

    in.exceptions(in.badbit | in.failbit | in.eofbit);

    kmers_.clear();
    bits_.clear();

    try {
        k_ = load_number(in);

        Deserializer deserializer(in);

        uint64_t tag = load_number(in);

        if (tag < std::numeric_limits<uint64_t>::max()) {
            packed_serialization_ = true;

            const auto size = tag;
            kmers_.reserve(size + 1);
            for (uint64_t i = 0; i < size; ++i) {
                kmers_.insert(deserializer.operator()<KmerPrefix>());
            }

            uint64_t size2 = load_number(in);
            values_.reserve(size2 + 1);
            for (uint64_t i = 0; i < size2; ++i) {
                values_.push_back(deserializer.operator()<std::array<uint64_t, 4>>());
            }

            //uint64_t size3 = load_number(in);
            //shifts_.reserve(size3 + 1);
            //for (uint64_t i = 0; i < size3; ++i) {
            //    shifts_.push_back(deserializer.operator()<bool>());
            //}

            uint64_t size3 = load_number(in);
            may_contain_source_kmer_.reserve(size3 + 1);
            for (uint64_t i = 0; i < size3; ++i) {
                may_contain_source_kmer_.push_back(deserializer.operator()<bool>());
            }

        } else {
            throw std::runtime_error("not implemented");
            //packed_serialization_ = false;
            //kmers_ = KmerIndex::deserialize(deserializer, true);
        }

        canonical_mode_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (!dynamic_cast<const DBGHashFast4Impl*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashFast4Impl*>(&other);

    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(canonical_mode_ == other_hash.canonical_mode_);
    assert(kmers_.size() == other_hash.kmers_.size());

    return kmers_ == other_hash.kmers_;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::get_index(const KmerConstIterator &kmer_iter) const {
    return ((kmer_iter - kmers_.begin()) << int(KMER::kBitsPerChar * 4)) + 1;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::get_index(const Kmer &kmer, const KmerConstIterator &kmer_iter) const {
    size_t shift = prefixes_[size_t(KmerPrefix(kmer.data()) & prefix_mask_)];

    KmerPrefix key = KmerPrefix(kmer.data()) & neither_mask_[shift];
    key = key >> int(KMER::kBitsPerChar * shift);
    KmerConstIterator iter = kmer_iter;

    if (iter == kmers_.end())
        return npos;

    const size_t rest = get_rest(kmer, shift);

    const uint64_t val = 1ull << rest;
    if (values_[iter - kmers_.begin()][shift] & val) {
        node_index base_node_index = ((iter - kmers_.begin()) << (KMER::kBitsPerChar * 4)) + 1;
        return base_node_index + (shift << (KMER::kBitsPerChar * 3)) + rest;
    }

    return npos;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::node_index
DBGHashFast4Impl<KMER>::get_index(const Kmer &kmer) const {
    size_t shift = prefixes_[size_t(KmerPrefix(kmer.data()) & prefix_mask_)];

    KmerPrefix key = KmerPrefix(kmer.data()) & neither_mask_[shift];
    key = key >> int(KMER::kBitsPerChar * shift);
    KmerConstIterator iter = kmers_.find(key);

    if (iter == kmers_.end())
        return npos;

    const size_t rest = get_rest(kmer, shift);

    const uint64_t val = 1ull << rest;
    if (values_[iter - kmers_.begin()][shift] & val) {
        node_index base_node_index = ((iter - kmers_.begin()) << (KMER::kBitsPerChar * 4)) + 1;
        return base_node_index + (shift << (KMER::kBitsPerChar * 3)) + rest;
    }

    return npos;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::KmerIterator
DBGHashFast4Impl<KMER>::get_iter(node_index node) {
    KmerIterator iter = kmers_.begin() + ((node - 1) >> (KMER::kBitsPerChar * 4));
    return iter;
}

template <typename KMER>
typename DBGHashFast4Impl<KMER>::KmerConstIterator
DBGHashFast4Impl<KMER>::get_const_iter(node_index node) const {
    KmerConstIterator iter = kmers_.begin() + ((node - 1) >> (KMER::kBitsPerChar * 4));
    return iter;
}

template <typename KMER>
const KMER DBGHashFast4Impl<KMER>
::get_kmer(node_index node,
           KmerConstIterator kmer_iter) const {
    assert(in_graph(node));

    KmerPrefix key = kmer_iter.key();
    KmerPrefix rest = KmerPrefix(node - 1) & first3_mask_;

    size_t shift = find_shift(node);

    KmerPrefix kmer = key << int(KMER::kBitsPerChar * shift);
    kmer |= rest >> int(KMER::kBitsPerChar * (3 - shift));
    kmer |= (rest & ~KmerPrefix(~0u << int(KMER::kBitsPerChar * (3 - shift)))) << int(KMER::kBitsPerChar * (k_ - (3 - shift)));
    return KMER(kmer);
}

template <typename KMER>
bool DBGHashFast4Impl<KMER>::in_graph(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    auto it = get_const_iter(node);
    KmerPrefix rest = KmerPrefix(node - 1) & first3_mask_;
    size_t shift = find_shift(node);
    return 1 & (values_[it - kmers_.begin()][shift] >> rest);
}

template <typename KMER>
void DBGHashFast4Impl<KMER>::call_nodes(const std::function<void(node_index)> &callback,
                               const std::function<bool()> &stop_early) const {
    for (auto iter = kmers_.begin(); iter != kmers_.end() && !stop_early(); ++iter) {
        for (size_t shift = 0; shift <= 3; ++shift) {
            auto val = values_[iter - kmers_.begin()][shift];
            while (val) {
                size_t i = __builtin_ffsl(val) - 1;
                if (val & (1ull << i)) {
                    auto node = get_index(iter) + (shift << (KMER::kBitsPerChar * 3)) + i;
                    assert(node <= num_nodes());
                    assert(node == kmer_to_node(get_node_sequence(node)));
                    callback(node);
                }
                val &= ~(1ull << i);
            }
        }
    }
}

std::unique_ptr<DBGHashFast4::DBGHashFast4Interface>
DBGHashFast4::initialize_graph(size_t k,
                                 bool canonical_mode,
                                 bool packed_serialization) {
    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashFast4Impl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode, packed_serialization
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashFast4Impl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode, packed_serialization
        );
    } else {
        return std::make_unique<DBGHashFast4Impl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode, packed_serialization
        );
    }
}

DBGHashFast4::DBGHashFast4(size_t k,
                               bool canonical_mode,
                               bool packed_serialization) {
    hash_dbg_ = initialize_graph(k, canonical_mode, packed_serialization);
}

bool DBGHashFast4::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        auto pos = in.tellg();
        auto k = load_number(in);
        in.seekg(pos, in.beg);

        // the actual value of |canonical| will be set in load
        hash_dbg_ = initialize_graph(k, false, false);
        return hash_dbg_->load(in) && in.good();
    } catch (...) {
        return false;
    }
}

bool DBGHashFast4::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}
