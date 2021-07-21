#include "boss.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>
#include <optional>
#include <stack>
#include <string>
#include <vector>

#include <progress_bar.hpp>
#include <tsl/hopscotch_set.h>

#include "common/threads/threading.hpp"
#include "common/serialization.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "boss_construct.hpp"


namespace mtg {
namespace graph {
namespace boss {

using common::logger;

#define CHECK_INDEX(idx) \
    assert(idx < W_->size()); \
    assert(idx > 0)
#define CHECK_NODE(idx) \
    assert(idx <= num_nodes()); \
    assert(idx > 0)

typedef BOSS::node_index node_index;
typedef BOSS::edge_index edge_index;
typedef BOSS::TAlphabet TAlphabet;

const size_t MAX_ITER_WAVELET_TREE_FAST = 1000;
const size_t MAX_ITER_WAVELET_TREE_DYN = 6;
const size_t MAX_ITER_WAVELET_TREE_STAT = 20;
const size_t MAX_ITER_WAVELET_TREE_SMALL = 1;

static const uint64_t kBlockSize = 9'999'872;
static_assert(!(kBlockSize & 0xFF));

const size_t TRAVERSAL_START_BATCH_SIZE = 100;
const size_t MAX_EDGE_QUEUE_SIZE = 10'000;
const size_t TASK_POOL_SIZE = 10'000;
const size_t MAX_DECODED_EDGES_IN_QUEUE = 20'000'000; // ~1 GB max


BOSS::BOSS(size_t k)
      : alph_size(kmer_extractor_.alphabet.size()),
        alphabet(kmer_extractor_.alphabet),
        bits_per_char_W_(sdsl::bits::hi(alph_size - 1) + 2),
        k_(k),
        last_(new bit_vector_dyn()),
        F_(alph_size, 0),
        NF_(alph_size, 0),
        W_(new wavelet_tree_dyn(bits_per_char_W_)) {

    assert(bits_per_char_W_ <= sizeof(TAlphabet) * 8
            && "Choose type for TAlphabet properly");

    assert(get_state() == BOSS::State::DYN);

    dynamic_cast<bit_vector_dyn&>(*last_).insert_bit(0, false);
    dynamic_cast<wavelet_tree_dyn&>(*W_).insert(0, 0);

    // add the dummy source node
    dynamic_cast<bit_vector_dyn&>(*last_).insert_bit(1, true);
    dynamic_cast<wavelet_tree_dyn&>(*W_).insert(0, 0);
    for (size_t j = 1; j < alph_size; j++) {
        F_[j] = 1;
        NF_[j] = 1;
    }
    assert(is_valid());
}

BOSS::BOSS(BOSSConstructor *builder) : BOSS::BOSS() {
    assert(builder);

    builder->build_graph(this);
    assert(is_valid());
}

BOSS::~BOSS() {
    delete W_;
    delete last_;
}

void BOSS::initialize(Chunk *chunk) {
    delete W_;
    delete last_;

    // TODO: optimize
    W_ = new wavelet_tree_stat(chunk->get_W_width(), chunk->W_);

    {
        chunk->last_.flush();
        sdsl::bit_vector last;
        std::ifstream in(chunk->last_.filename(), std::ios::binary);
        last.load(in);
        last_ = new bit_vector_stat(std::move(last));
    }

    F_ = chunk->F_;
    recompute_NF();

    k_ = chunk->k_;
    // TODO:
    // alph_size = chunk->alph_size_;

    state = State::STAT;
}

/**
 * Given a pointer to BOSS tables G1 and G2, the function compares their elements to the
 * each other. It will perform an element wise comparison of the arrays W, last and
 * F and will only check for identity. If any element differs, the function will return
 * false and true otherwise.
 */
bool BOSS::equals_internally(const BOSS &other, bool verbose) const {
    // compare size
    if (num_edges() != other.num_edges()) {
        if (verbose)
            std::cout << "sizes of graphs differ"
                      << "\n1: " << W_->size()
                      << "\n2: " << other.W_->size()
                      << std::endl;
        return false;
    }

    assert(F_.size() == other.F_.size());

    bool all_equal = (*W_ == *other.W_
                        && *last_ == *other.last_
                        && F_ == other.F_);

    if (all_equal || !verbose)
        return all_equal;

    // compare last
    for (edge_index i = 0; i < W_->size(); ++i) {
        if (get_last(i) != other.get_last(i)) {
            if (verbose)
                std::cout << "last differs at position " << i
                          << "\n1: last[" << i << "] = " << get_last(i)
                          << "\n2: last[" << i << "] = " << other.get_last(i)
                          << std::endl;
            return false;
        }
    }

    // compare W
    for (edge_index i = 0; i < W_->size(); ++i) {
        if (get_W(i) != other.get_W(i)) {
            if (verbose)
                std::cout << "W differs at position " << i
                          << "\n1: W[" << i << "] = " << get_W(i)
                          << "\n2: W[" << i << "] = " << other.get_W(i)
                          << std::endl;
            return false;
        }
    }

    // compare F
    for (size_t c = 0; c < F_.size(); ++c) {
        if (get_F(c) != other.get_F(c)) {
            if (verbose)
                std::cout << "F differs at position " << c
                          << "\n1: F[" << c << "] = " << get_F(c)
                          << "\n2: F[" << c << "] = " << other.get_F(c)
                          << std::endl;
            return false;
        }
    }

    assert(all_equal);
    return true;
}

/**
 * Check whether BOSS tables store the same data.
 * FYI: this function reconstructs all the kmers, so
 * the complexity is at least O(k x n).
 */
bool BOSS::operator==(const BOSS &other) const {
    edge_index i = 1;
    edge_index j = 1;

    while (i < W_->size() && j < other.W_->size()) {

        std::string first_node, second_node;
        bool first_last, second_last;
        char first_label, second_label;

        do {
            first_node = get_node_str(i);
            first_last = get_last(i);
            first_label = decode(get_W(i) % alph_size);
            i++;
        } while (first_node.find(kSentinel) != std::string::npos && i < W_->size());

        do {
            second_node = other.get_node_str(j);
            second_last = other.get_last(j);
            second_label = decode(other.get_W(j) % alph_size);
            j++;
        } while (second_node.find(kSentinel) != std::string::npos && j < other.W_->size());

        if (i == W_->size() || j == other.W_->size())
            break;

        if (first_node != second_node
                || first_last != second_last
                || first_label != second_label)
            return false;
    }
    return i == W_->size() && j == other.W_->size();
}

void BOSS::serialize(const std::string &filename) const {
    const auto out_filename = utils::make_suffix(filename, kExtension);

    std::ofstream outstream(out_filename, std::ios::binary);
    if (!outstream.good()) {
        throw std::ofstream::failure(
            std::string("Error: Can't write to file ") + out_filename
        );
    }

    serialize(outstream);

    outstream.close();
}

void BOSS::serialize(std::ofstream &outstream) const {
    if (!outstream.good())
        throw std::ofstream::failure("Error: Can't write to file");

    // write F values, k, and state
    serialize_number_vector_raw(outstream, F_);
    serialize_number(outstream, k_);
    serialize_number(outstream, state);
    outstream.flush();

    // write Wavelet Tree
    W_->serialize(outstream);
    outstream.flush();

    // write last array
    last_->serialize(outstream);
    outstream.flush();
}

bool BOSS::load(const std::string &filename) {
    auto file = utils::make_suffix(filename, kExtension);

    std::ifstream instream(file, std::ios::binary);

    return load(instream);
}

bool BOSS::load(std::ifstream &instream) {
    // if not specified in the file, the default for loading is dynamic
    state = State::DYN;

    try {
        // load F, k, and state
        F_ = load_number_vector_raw<edge_index>(instream);
        k_ = load_number(instream);
        state = static_cast<State>(load_number(instream));

        if (F_.size() != alph_size) {
            std::cerr << "ERROR: failed to load F vector, incompatible size" << std::endl;
            return false;
        }

        // load W and last arrays
        delete W_;
        delete last_;
        switch (state) {
            case State::DYN:
                W_ = new wavelet_tree_dyn(bits_per_char_W_);
                last_ = new bit_vector_dyn();
                break;
            case State::STAT:
                W_ = new wavelet_tree_stat(bits_per_char_W_);
                last_ = new bit_vector_stat();
                break;
            case State::FAST:
                W_ = new wavelet_tree_fast(bits_per_char_W_);
                last_ = new bit_vector_stat();
                break;
            case State::SMALL:
                W_ = new wavelet_tree_small(bits_per_char_W_);
                last_ = new bit_vector_small();
                break;
        }
        if (!W_->load(instream)) {
            std::cerr << "ERROR: failed to load W vector" << std::endl;
            return false;
        }

        if (!last_->load(instream)) {
            std::cerr << "ERROR: failed to load L vector" << std::endl;
            return false;
        }

        recompute_NF();

        return instream.good();

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load the BOSS table." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void BOSS::serialize_suffix_ranges(std::ofstream &outstream) const {
    // dump node range index
    serialize_number(outstream, indexed_suffix_length_);

    outstream.write(reinterpret_cast<const char *>(indexed_suffix_ranges_.data()),
                    indexed_suffix_ranges_.size()
                        * sizeof(decltype(indexed_suffix_ranges_)::value_type));
}

bool BOSS::load_suffix_ranges(std::ifstream &instream) {
    // load node suffix range index if exists
    try {
        indexed_suffix_length_ = load_number(instream);
        if (!indexed_suffix_length_
                || indexed_suffix_length_ > k_
                || indexed_suffix_length_ * log2(alph_size - 1) > 63)
            throw std::ifstream::failure("");

        uint64_t index_size = 1;
        for (size_t len = 1; len <= indexed_suffix_length_; ++len) {
            index_size *= (alph_size - 1);
        }

        indexed_suffix_ranges_.resize(index_size);

        instream.read(reinterpret_cast<char *>(indexed_suffix_ranges_.data()),
                      indexed_suffix_ranges_.size()
                        * sizeof(decltype(indexed_suffix_ranges_)::value_type));

        if (!instream.good())
            throw std::ifstream::failure("");

        return true;

    } catch(...) {
        indexed_suffix_length_ = 0;
        indexed_suffix_ranges_.clear();
        return false;
    }
}

//
//
// HELPER QUERY FUNCTIONS
//
//

/**
 * For the given position i in W and a character c from the alphabet,
 * return the number of occurrences of c in W up to (including) position i.
 */
uint64_t BOSS::rank_W(edge_index i, TAlphabet c) const {
    assert(i < W_->size());

    return i == 0 ? 0 : W_->rank(c, i) - (c == 0);
}

/**
 * Return the position of the last occurrence of |c| in W[1..i] or zero
 * if such does not exist.
 */
edge_index BOSS::pred_W(edge_index i, TAlphabet c) const {
    CHECK_INDEX(i);

    edge_index prev = W_->prev(i, c);
    return prev < W_->size() ? prev : 0;
}

/**
 * For characters |first| and |second|, return the last occurrence
 * of them in W[1..i], i.e. max(pred_W(i, first), pred_W(i, second)).
 */
edge_index BOSS::pred_W(edge_index i, TAlphabet c_first, TAlphabet c_second) const {
    CHECK_INDEX(i);
    assert(c_first != c_second);

    // trying to avoid calls of W_->prev
    uint64_t max_iter = 0;
    switch (state) {
        case STAT:
            max_iter = MAX_ITER_WAVELET_TREE_STAT;
            break;
        case DYN:
            max_iter = MAX_ITER_WAVELET_TREE_DYN;
            break;
        case SMALL:
            max_iter = MAX_ITER_WAVELET_TREE_SMALL;
            break;
        case FAST:
            max_iter = MAX_ITER_WAVELET_TREE_FAST;
            break;
    }

    edge_index end = i - std::min(i, max_iter);

    // TODO: add a cap -- iterate alph_size^2 elements at maximum
    while (i > end) {
        TAlphabet w = get_W(i);
        if (w == c_first || w == c_second)
            return i;
        i--;
    }

    if (!i)
        return 0;

    // get the previous position via rank + select calls
    uint64_t r_first = W_->rank(c_first, i);
    edge_index select_first = r_first ? W_->select(c_first, r_first) : 0;

    uint64_t r_second = W_->rank(c_second, i);
    edge_index select_second = r_second ? W_->select(c_second, r_second) : 0;

    return std::max(select_first, select_second);
}

/**
 * Return the position of the first occurrence of |c| in W[i..N].
 */
edge_index BOSS::succ_W(edge_index i, TAlphabet c) const {
    CHECK_INDEX(i);

    return W_->next(i, c);
}

/**
 * For characters |first| and |second|, return the first occurrence
 * of them in W[i..N], i.e. min(succ_W(i, first), succ_W(i, second)).
 */
std::pair<edge_index, TAlphabet>
BOSS::succ_W(edge_index i, TAlphabet c_first, TAlphabet c_second) const {
    CHECK_INDEX(i);
    assert(c_first != c_second);

    // trying to avoid calls of W_->next
    uint64_t max_iter = 0;
    switch (state) {
        case STAT:
            max_iter = MAX_ITER_WAVELET_TREE_STAT;
            break;
        case DYN:
            max_iter = MAX_ITER_WAVELET_TREE_DYN;
            break;
        case SMALL:
            max_iter = MAX_ITER_WAVELET_TREE_SMALL;
            break;
        case FAST:
            max_iter = MAX_ITER_WAVELET_TREE_FAST;
            break;
    }

    edge_index end = i + std::min(W_->size() - i, max_iter);

    while (i < end) {
        TAlphabet w = get_W(i);
        if (w == c_first)
            return std::make_pair(i, c_first);
        if (w == c_second)
            return std::make_pair(i, c_second);
        i++;
    }

    if (i == W_->size())
        return std::make_pair(W_->size(), 0);

    // get the next position via rank + select calls
    uint64_t r_first = W_->rank(c_first, i - 1) + 1;
    edge_index select_first = r_first <= W_->count(c_first)
                                ? W_->select(c_first, r_first)
                                : W_->size();
    uint64_t r_second = W_->rank(c_second, i - 1) + 1;
    edge_index select_second = r_second <= W_->count(c_second)
                                ? W_->select(c_second, r_second)
                                : W_->size();

    if (select_first < select_second) {
        return std::make_pair(select_first, c_first);
    } else if (select_second < select_first) {
        return std::make_pair(select_second, c_second);
    } else {
        assert(select_first == W_->size());
        assert(select_second == W_->size());
        return std::make_pair(W_->size(), 0);
    }
}

/**
 * Transforms a boss edge index to the index of its source node.
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
node_index BOSS::rank_last(edge_index i) const {
    assert(i < last_->size());

    return i == 0 ? 0 : last_->rank1(i);
}

/**
 * Transforms a boss node index to the index of its last outgoing edge.
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
edge_index BOSS::select_last(node_index i) const {
    assert(i <= last_->num_set_bits());

    return i == 0 ? 0 : last_->select1(i);
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
edge_index BOSS::pred_last(edge_index i) const {
    if (!i)
        return 0;

    CHECK_INDEX(i);

    edge_index prev = last_->prev1(i);

    return prev < last_->size() ? prev : 0;
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
edge_index BOSS::succ_last(edge_index i) const {
    CHECK_INDEX(i);

    return last_->next1(i);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character.
 */
edge_index BOSS::bwd(edge_index i) const {
    CHECK_INDEX(i);

    node_index target_node = rank_last(i - 1) + 1;

    if (target_node == 1)
        return 1;

    // get value of last position in node i
    TAlphabet c = get_node_last_value(i);
    assert(c && "there must be no edges of type ***$* except for the main dummy");
    // get the offset among the representative edges ending with |c| and select
    return W_->select(c, target_node - NF_[c]);
}

/**
 * This function gets a position i reflecting the r-th occurrence of the corresponding
 * character c in W and returns the position of the r-th occurrence of c in last.
 */
edge_index BOSS::fwd(edge_index i, TAlphabet c) const {
    CHECK_INDEX(i);

    // |c| must be the value in W at position i
    assert(c == get_W(i) % alph_size);
    assert(i == 1 || c != kSentinelCode);
    // get the index of the target node
    node_index target_node = NF_[c] + rank_W(i, c);
    // get the index of the representative edge with that target node
    return select_last(target_node);
}

/**
 * Using the offset structure F this function returns the value of the last
 * position of the source node for edge i.
 */
TAlphabet BOSS::get_node_last_value(edge_index i) const {
    CHECK_INDEX(i);

    if (i == 0)
        return 0;

    for (TAlphabet c = 0; c < alph_size; c++) {
        if (F_[c] >= i)
            return c - 1;
    }
    return alph_size - 1;
}

/**
 * Given index i of an edge and a value k, this function
 * returns the k-th last character of the source node for edge i.
 */
std::pair<TAlphabet, edge_index>
BOSS::get_minus_k_value(edge_index i, size_t k) const {
    CHECK_INDEX(i);

    for (; k > 0; --k) {
        i = bwd(i);
    }
    return std::make_pair(get_node_last_value(i), i);
}

/**
 * Given an edge index |i| and a character |c|, get the index of the edge with
 * label c outgoing from the same source node if such exists and npos otherwise.
 */
edge_index BOSS::pick_edge(edge_index edge, TAlphabet c) const {
    CHECK_INDEX(edge);
    assert(get_last(edge) && "must be the last outgoing edge");
    assert(c < alph_size);

    do {
        TAlphabet w = get_W(edge);
        if (w == c || w == c + alph_size)
            return edge;
    } while (--edge && !get_last(edge));

    return npos;
}

/**
 * Given an edge index |x| and a character |c|, get the index of an adjacent
 * incoming edge with the first character c if such exists and npos otherwise.
 */
edge_index BOSS::pick_incoming_edge(edge_index x, TAlphabet c) const {
    CHECK_INDEX(x);
    assert(get_W(x) < alph_size && "must be the first incoming edge");
    assert(c <= alph_size);

    if (c == alph_size)
        return npos;

    // only one incoming edge for the dummy source node
    if (x == 1) {
        if (c == kSentinelCode) {
            return 1;
        } else {
            return npos;
        }
    }

    // check if the first incoming edge has label `c`
    if (get_minus_k_value(x, k_ - 1).first == c)
        return x;

    if (x + 1 == W_->size())
        return npos;

    // TODO: could be improved. implement without succ_W.
    TAlphabet d = get_W(x);
    edge_index y = succ_W(x + 1, d);

    // iterate over the rest of the incoming edges
    while (x + 1 < y) {
        x = succ_W(x + 1, d + alph_size);
        if (x < y && get_minus_k_value(x, k_ - 1).first == c) {
            return x;
        }
    }
    return npos;
}

void BOSS::call_incoming_to_target(edge_index edge, TAlphabet d,
                                   const Call<edge_index> &callback) const {
    CHECK_INDEX(edge);
    assert(get_W(edge) < alph_size && "must be the first incoming edge");
    assert(d == get_W(edge));

    callback(edge);

    // iterate through all indices with edge label d + alph_size
    // which are less than the next index with edge label d
    TAlphabet d_next;
    while (++edge < W_->size()) {

        std::tie(edge, d_next) = succ_W(edge, d, d + alph_size);

        if (d_next != d + alph_size)
            break;

        callback(edge);
    }
}

/**
 * Given an edge index i, this function returns true if that is
 * the only outgoing edge from its source node.
 */
bool BOSS::is_single_outgoing(edge_index i) const {
    CHECK_INDEX(i);

    return get_last(i) && (i == 1 || get_last(i - 1));
}

/**
 * Given an edge index i and label w, this function returns
 * true if that is the only edge incoming to its target node.
 */
bool BOSS::is_single_incoming(edge_index i, TAlphabet w) const {
    CHECK_INDEX(i);
    assert(w == get_W(i));
    assert(w != alph_size);

    if (w > alph_size)
        return false;

    // start from the next edge
    i++;

    return i == W_->size()
            || succ_W(i, w, w + alph_size).second != w + alph_size;
}

/**
 * Given an edge index i (first incoming) and label d, this function returns
 * the number of edges incoming to its target node.
 */
size_t BOSS::num_incoming_to_target(edge_index x, TAlphabet d) const {
    CHECK_INDEX(x);
    assert(get_W(x) < alph_size && "must be the first incoming edge");
    assert(d == get_W(x));

    if (x + 1 == W_->size())
        return 1;

#if 0 // the second section is faster for all existing graph states
    edge_index y = succ_W(x + 1, d);
    return 1 + rank_W(y - 1, d + alph_size) - rank_W(x - 1, d + alph_size);
#else
    size_t indeg = 0;
    call_incoming_to_target(x, d, [&indeg](auto) { indeg++; });
    assert(indeg && "there is always at least one incoming edge");
    return indeg;
#endif
}


/**
 * Given a node label kmer, this function returns the index
 * of the corresponding node or the closest predecessor, if no node
 * with the sequence is not found.
 */
node_index BOSS::pred_kmer(const std::vector<TAlphabet> &kmer) const {
    assert(kmer.size() == k_);

    // get first
    auto kmer_it = kmer.begin();

    edge_index last_ = *kmer_it + 1 < alph_size
                        ? F_[*kmer_it + 1]
                        : W_->size() - 1;
    uint64_t shift = 0;

    // update range iteratively while scanning through s
    while (++kmer_it != kmer.end()) {
        TAlphabet s = *kmer_it;
        assert(s <= alph_size);

        // invalid character that does not belong to the alphabet
        // assume |s| greater than all characters
        if (s == alph_size) {
            last_ = W_->size() - 1;
            shift = 0;
            continue;
        }

        edge_index last_target = pred_W(last_, s, s + alph_size);
        if (last_target > 0) {
            if (rank_last(last_target - 1) < rank_last(last_ - 1))
                shift = 0;
            last_ = fwd(last_target, s);
            continue;
        }
        assert(s > 0);

        last_target = succ_W(last_, s, s + alph_size).first;

        if (last_target < W_->size()) {
            last_ = fwd(last_target, s);
            shift = 1;
        } else {
            last_ = F_[s];
            shift = 0;
        }
    }

    CHECK_NODE(rank_last(last_ - shift));
    return rank_last(last_ - shift);
}


/**
 * This function gets two edge indices and returns if their source
 * node labels share a k-1 suffix.
 */
bool BOSS::compare_node_suffix(edge_index first, edge_index second) const {
    for (size_t i = 0; i < k_ - 1; ++i) {
        if (first == second)
            return true;

        if (get_node_last_value(first) != get_node_last_value(second)) {
            return false;
        }
        first = bwd(first);
        second = bwd(second);
    }
    return true;
}

/**
 * This function gets an edge indix and checks if its source
 * node has the same k-1 suffix as k-mer |second|.
 */
bool BOSS::compare_node_suffix(edge_index first, const TAlphabet *second) const {
    for (auto it = second + k_ - 1; it > second; --it) {
        if (get_node_last_value(first) != *it) {
            return false;
        }
        first = bwd(first);
    }
    return true;
}

/**
 * Given an edge index i, this function returns the k-mer sequence of its
 * source node, and the node whose last character corresponds to the first
 * character of the sequence. If the graph is suffix indexed, then the returned
 * node is the last node visited after k - indexed_suffix_length_ bwd steps.
 */
std::tuple<std::vector<TAlphabet>, edge_index, bool> BOSS
::get_node_seq_with_end_node_indexed(edge_index x) const {
    CHECK_INDEX(x);

    std::vector<TAlphabet> ret(k_);
    size_t i = k_;

    if (indexed_suffix_length_) {
        while (i > indexed_suffix_length_) {
            CHECK_INDEX(x);

            ret[--i] = get_node_last_value(x);
            x = bwd(x);
        }

        auto it = std::lower_bound(
            indexed_suffix_ranges_.begin(),
            indexed_suffix_ranges_.end(),
            x,
            [](const auto &range, edge_index edge) { return (range.second < edge); }
        );

        if (it != indexed_suffix_ranges_.end() && it->first <= x) {
            assert(x <= it->second);
            uint64_t index = it - indexed_suffix_ranges_.begin();
            for (i = 0; i < indexed_suffix_length_; ++i) {
                uint64_t next_index = index / (alph_size - 1);
                ret[i] = index - next_index * (alph_size - 1) + 1;
                index = next_index;
            }
            return std::make_tuple(ret, x, true);
        }
    }

    ret[--i] = get_node_last_value(x);

    while (i > 0) {
        CHECK_INDEX(x);

        x = bwd(x);
        ret[--i] = get_node_last_value(x);
    }

    return std::make_tuple(ret, x, false);
}

/**
 * Given an edge index i, this function returns the k-mer sequence of its
 * source node.
 */
std::vector<TAlphabet> BOSS::get_node_seq(edge_index i) const {
    return std::get<0>(get_node_seq_with_end_node_indexed(i));
}

/**
 * Given an edge index i, this function returns the k-mer sequence of its
 * source node, and the node whose last character corresponds to the first
 * character of the sequence.
 */
std::pair<std::vector<TAlphabet>, edge_index> BOSS
::get_node_seq_with_end_node(edge_index i) const {
    auto [seq, last_node, indexed] = get_node_seq_with_end_node_indexed(i);

    if (indexed) {
        assert(indexed_suffix_length_);

        auto it = seq.rend() - indexed_suffix_length_;
        assert(*it == get_node_last_value(last_node));
        for (++it; it != seq.rend(); ++it) {
            last_node = bwd(last_node);
            assert(*it == get_node_last_value(last_node));
        }
    }

    assert(seq[0] == get_node_last_value(last_node));
    assert(std::make_pair(seq[0], last_node) == get_minus_k_value(i, k_ - 1));

    return std::make_pair(std::move(seq), last_node);
}

/**
 * Given a node index k_node, this function returns the k-mer sequence of the
 * node as a string.
 */
std::string BOSS::get_node_str(edge_index k_node) const {
    CHECK_INDEX(k_node);
    return decode(get_node_seq(k_node));
}

void BOSS::map_to_edges(std::string_view sequence,
                        const std::function<void(edge_index)> &callback,
                        const std::function<bool()> &terminate,
                        const std::function<bool()> &skip) const {
    map_to_edges(encode(sequence), callback, terminate, skip);
}

std::vector<edge_index> BOSS::map_to_edges(std::string_view sequence) const {
    return map_to_edges(encode(sequence));
}

void BOSS::map_to_edges(const std::vector<TAlphabet> &seq_encoded,
                        const std::function<void(edge_index)> &callback,
                        const std::function<bool()> &terminate,
                        const std::function<bool()> &skip) const {
    assert(std::all_of(seq_encoded.begin(), seq_encoded.end(),
                       [this](TAlphabet c) { return c <= alph_size; }));

    if (seq_encoded.size() <= k_)
        return;

    // Mark where (k+1)-mers with invalid characters end
    // Example for (k+1)=3: [X]***[X]****[X]***
    //              ---->   [111]0[111]00[111]0
    auto invalid = utils::drag_and_mark_segments(seq_encoded, alph_size, k_ + 1);

    // slide through all (k+1)-mers
    for (size_t i = 0; i + k_ + 1 <= seq_encoded.size() && !terminate(); ++i) {
        if (skip())
            continue;

        if (invalid[i + k_]) {
            // this (k+1)-mer contains at least one invalid character
            callback(npos);
            continue;
        }

        auto edge = map_to_edge(seq_encoded.data() + i,
                                seq_encoded.data() + i + k_ + 1);
        callback(edge);

        while (edge && ++i + k_ < seq_encoded.size()) {
            if (terminate())
                return;

            if (skip())
                break;

            if (invalid[i + k_]) {
                // this (k+1)-mer contains at least one invalid character
                callback(npos);
                break;
            }

            edge = fwd(edge, seq_encoded[i + k_ - 1]);
            edge = pick_edge(edge, seq_encoded[i + k_]);

            callback(edge);
        }
    }
}

std::vector<edge_index>
BOSS::map_to_edges(const std::vector<TAlphabet> &seq_encoded) const {
    std::vector<edge_index> indices;
    indices.reserve(seq_encoded.size());
    map_to_edges(seq_encoded,
                 [&indices](edge_index i) { indices.push_back(i); });
    return indices;
}

/**
 * Returns the number of nodes in BOSS graph.
 */
uint64_t BOSS::num_nodes() const {
    return last_->num_set_bits();
}

/**
 * Return the number of edges in BOSS graph.
 */
uint64_t BOSS::num_edges() const {
    return W_->size() - 1;
}

/**
 * This function gets a character |c| and updates the edge offsets in F_
 * by incrementing them with +1 (for edge insertion) or decrementing
 * (for edge delition) and updates the node offsets NF_ accordingly.
 * The node offsets are updated in two cases:
 *  - the edge was inserted (delta = +1) and it is the only edge outgoing
 *  from its source node,
 *  - the edge was erased (delta = -1) and its source node has no other
 *  outgoing edges.
 * Pass |is_representative| = `true` if one of these two conditions is
 * satisfied and `false` otherwise.
 */
void BOSS::update_F(TAlphabet c, int delta, bool is_representative) {
    assert(c < alph_size);
    assert(std::abs(delta) == 1);

    for (TAlphabet i = c + 1; i < alph_size; i++) {
        F_[i] += delta;
        NF_[i] += delta * is_representative;
    }
}

/**
 * Recompute the node offsets NF_ from F_. Call after changes in F_.
 */
void BOSS::recompute_NF() {
    // precompute the node offsets
    NF_.resize(F_.size());
    for (size_t c = 0; c < F_.size(); ++c) {
        NF_[c] = rank_last(F_[c]);
    }
}

TAlphabet BOSS::encode(char s) const {
    assert(kmer_extractor_.encode(kSentinel) != kSentinelCode);
    assert(kmer_extractor_.encode(s) <= alph_size);
    return kmer_extractor_.encode(s);
}

std::vector<TAlphabet> BOSS::encode(std::string_view sequence) const {
    std::vector<TAlphabet> seq_encoded = kmer_extractor_.encode(sequence);
    assert(std::all_of(seq_encoded.begin(), seq_encoded.end(),
                       [this](TAlphabet c) { return c <= alph_size; }));
    return seq_encoded;
}

char BOSS::decode(TAlphabet c) const {
    assert(kmer_extractor_.alphabet[kSentinelCode] == kSentinel);
    assert(c < alph_size);
    return kmer_extractor_.decode(c);
}

std::string BOSS::decode(const std::vector<TAlphabet> &seq_encoded) const {
    assert(kmer_extractor_.encode(kSentinel) != kSentinelCode);
    assert(std::all_of(seq_encoded.begin(), seq_encoded.end(),
                       [this](TAlphabet c) { return c < alph_size; }));
    return kmer_extractor_.decode(seq_encoded);
}

template <class WaveletTree, class BitVector>
void convert(wavelet_tree **W_, bit_vector **last_) {
    wavelet_tree *W_new = new WaveletTree((*W_)->convert_to<WaveletTree>());
    delete *W_;
    *W_ = W_new;

    bit_vector *last_new = new BitVector((*last_)->convert_to<BitVector>());
    delete *last_;
    *last_ = last_new;
}

void BOSS::switch_state(State new_state) {

    //std::cerr << "switching state from " << this->state << " to " << state << std::endl;
    if (state == new_state)
        return;

    switch (new_state) {
        case State::STAT: {
            convert<wavelet_tree_stat, bit_vector_stat>(&W_, &last_);
            break;
        }
        case State::SMALL: {
            convert<wavelet_tree_small, bit_vector_small>(&W_, &last_);
            break;
        }
        case State::FAST: {
            convert<wavelet_tree_fast, bit_vector_stat>(&W_, &last_);
            break;
        }
        case State::DYN: {
            convert<wavelet_tree_dyn, bit_vector_dyn>(&W_, &last_);
            break;
        }
    }
    state = new_state;
}

void BOSS::print_internal_representation(std::ostream &os) const {
    os << "F:";
    for (auto i : F_) {
        os << " " << i;
    }
    os << "\nL: " << *last_;
    os << "\nW:";
    for (edge_index i = 0; i < W_->size(); ++i) {
        os << " " << static_cast<int>(get_W(i));
    }
    os << std::endl;
}

void BOSS::print(std::ostream &os) const {
    assert(is_valid());

    std::string vertex_header = "Vertex";
    vertex_header.resize(k_, ' ');

    os << "Index" << "\t" << "L"
                  << "\t" << vertex_header
                  << "\t" << "W" << std::endl;

    for (edge_index i = 1; i < W_->size(); i++) {
        TAlphabet w = get_W(i);
        assert(w != alph_size);
        os << i << "\t" << get_last(i)
                << "\t" << get_node_str(i)
                << "\t" << decode(w % alph_size)
                        << (w > alph_size ? "-" : "")
                        << std::endl;
    }
}

void BOSS::print_adj_list(std::ostream &os) const {
    for (edge_index edge = 1; edge < W_->size(); ++edge) {
        os << 1 + rank_last(fwd(edge, get_W(edge) % alph_size) - 1)
           << " ";
        if (get_last(edge))
            os << "\n";
    }
}

///////////////
// Construct //
///////////////

// add a full sequence to the graph
void BOSS::add_sequence(std::string_view seq,
                        bool try_extend,
                        std::vector<edge_index> *edges_inserted) {
    if (seq.size() < k_ + 1)
        return;

    if (get_state() != State::DYN)
        throw std::runtime_error("representation must be dynamic");

    // prepend k buffer characters, in case we need to start with dummy node
    std::vector<TAlphabet> sequence(seq.size() + k_);
    std::transform(seq.begin(), seq.end(), sequence.begin() + k_,
        [this](char c) { return this->encode(c); }
    );

    TAlphabet *begin_segm = sequence.data() + k_;
    TAlphabet *end_segm;
    TAlphabet *end = sequence.data() + sequence.size();

    while (begin_segm + k_ < end) {

        assert(std::all_of(begin_segm, end,
                           [this](TAlphabet c) { return c <= alph_size; }));

        end_segm = std::find(begin_segm, end, alph_size);

        if (begin_segm + k_ < end_segm) {

            edge_index source;

            if (!try_extend || !(source = index(begin_segm, begin_segm + k_))) {
                // start insertion from the main dummy source node
                begin_segm -= k_;
                std::fill(begin_segm, begin_segm + k_, kSentinelCode);
                source = 1; // the dummy source node
            }

            while (begin_segm + k_ < end_segm) {
                source = append_pos(begin_segm[k_], source, begin_segm, edges_inserted);
                begin_segm++;
            }

            verbose_cout("edges ", num_edges(), " / nodes ", num_nodes(), "\n");
        }

        begin_segm = end_segm + 1;
    }
}

/**
 * Given a character c and an edge index, this function
 * creates an outgoing edge from the same source node with
 * label c if it is not a part of the graph yet.
 */
edge_index BOSS::append_pos(TAlphabet c, edge_index source_node,
                            const TAlphabet *source_node_kmer,
                            std::vector<edge_index> *edges_inserted) {
    CHECK_INDEX(source_node);
    assert(source_node_kmer);
    assert(std::vector<TAlphabet>(source_node_kmer, source_node_kmer + k_)
                                                == get_node_seq(source_node));
    assert(c < alph_size);
    assert(get_state() == State::DYN);

    // get range of identical nodes (without W) pos current end position
    edge_index begin = pred_last(source_node - 1) + 1;
    edge_index end = succ_last(source_node) + 1;

    // get position of the first occurrence of c or c- in W after p
    edge_index prev_c_pos = pred_W(end - 1, c, c + alph_size);
    // if the edge already exists, traverse it and return the index
    if (prev_c_pos >= begin)
        return fwd(prev_c_pos, c);

    /**
     * We found that c does not exist in the current range yet and now have to
     * figure out if we need to add c or c- to the range.
     * To do this, we check if there is a preceding position j1 with W[j1] == c
     * whose node shares a k-1 suffix with the current node.
     * If yes, we add c- instead of c.
     */

    // Check if the new edge will be the first incoming for its target node
    if (prev_c_pos > 0 && compare_node_suffix(prev_c_pos, source_node_kmer)) {
        // the new edge will not be the first incoming for its target node
        // insert the edge
        edge_index inserted = insert_edge(c + alph_size, begin, end);
        if (edges_inserted && inserted)
            edges_inserted->push_back(inserted);

        return fwd(prev_c_pos, c);
    }

    // The new edge will be the first incoming for its target node,
    // and therefore the new edge will be marked by c (not c-)

    // The inserted node may share its k-1 suffix with one of the next nodes
    // Get the position of the first occurrence of c after p (including p + 1)
    edge_index first_c = end < W_->size()
                            ? succ_W(end, c)
                            : W_->size();

    bool the_only_incoming = true;
    if (first_c < W_->size()
            && !(the_only_incoming = !compare_node_suffix(first_c, source_node_kmer))) {
        // The inserted edge will not be the only incoming for its target node.
        // Relabel the next incoming edge from c to c- since
        // the new edge is the first incoming for that target node.
        dynamic_cast<wavelet_tree_dyn&>(*W_).set(first_c, c + alph_size);
    }

    // insert the edge
    edge_index inserted = insert_edge(c, begin, end);
    if (edges_inserted && inserted)
        edges_inserted->push_back(inserted);

    if (!the_only_incoming)
        return fwd(first_c + (inserted > 0), c);

    // The inserted node forms a dead-end, thus a sentinel must be added
    edge_index sentinel_pos = select_last(NF_[c] + rank_W(begin - 1, c)) + 1;

    update_F(c, +1, true);
    dynamic_cast<wavelet_tree_dyn&>(*W_).insert(sentinel_pos, kSentinelCode);
    dynamic_cast<bit_vector_dyn&>(*last_).insert_bit(sentinel_pos, true);

    if (edges_inserted)
        edges_inserted->push_back(sentinel_pos);

    assert((*W_)[0] == 0);

    return sentinel_pos;
}


edge_index BOSS::insert_edge(TAlphabet c, edge_index begin, edge_index end) {
    assert(c != alph_size);
    assert(c < 2 * alph_size);
    assert(get_state() == State::DYN);

    if (begin > 1 && get_W(begin) == kSentinelCode) {
        // the source node is the dead-end with outgoing sentinel
        // replace this sentinel with the proper label
        dynamic_cast<wavelet_tree_dyn&>(*W_).set(begin, c);
        return 0;
    } else {
        // the source node already has some outgoing edges

        // find the exact position of the new edge
        edge_index pos = begin;
        while (pos < end && get_W(pos) % alph_size < c % alph_size) {
            pos++;
        }

        // insert the new edge
        update_F(get_node_last_value(begin), +1, false);
        dynamic_cast<bit_vector_dyn&>(*last_).insert_bit(begin, false);
        dynamic_cast<wavelet_tree_dyn&>(*W_).insert(pos, c);

        assert(pos);
        return pos;
    }
}


// Given an edge list, remove them from the BOSS graph.
// TODO: fix the implementation (anchoring the isolated nodes)
void BOSS::erase_edges_dyn(const std::set<edge_index> &edges) {
    uint64_t shift = 0;

    if (get_state() != State::DYN)
        throw std::runtime_error("representation must be dynamic");

    for (edge_index edge : edges) {
        assert(edge >= shift);
        edge_index edge_id = edge - shift;

        TAlphabet d = get_W(edge_id);
        if (d < alph_size && edge_id + 1 < W_->size()) {
            //fix W array
            auto [next, d_next] = succ_W(edge_id + 1, d, d + alph_size);
            if (d_next == d + alph_size)
                dynamic_cast<wavelet_tree_dyn&>(*W_).set(next, d);
        }
        dynamic_cast<wavelet_tree_dyn&>(*W_).remove(edge_id);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        bool last = get_last(edge_id);
        if (last && (edge >= shift + 1) && !get_last(edge_id - 1)) {
            update_F(get_node_last_value(edge_id), -1, false);
            dynamic_cast<bit_vector_dyn&>(*last_).delete_bit(edge_id - 1);
        } else {
            update_F(get_node_last_value(edge_id), -1, last);
            dynamic_cast<bit_vector_dyn&>(*last_).delete_bit(edge_id);
        }
        shift++;
    }
}

/**
 * Erase exactly all the masked edges from the BOSS graph,
 * may invalidate the BOSS table (if leaves nodes with no incoming edges).
 * Returns the number of edges erased.
 */
uint64_t BOSS::erase_edges(const sdsl::bit_vector &edges_to_remove_mask) {
    uint64_t num_edges_to_remove = sdsl::util::cnt_one_bits(edges_to_remove_mask);
    if (!num_edges_to_remove)
        return 0;

    // update last
    sdsl::bit_vector old_last = last_->convert_to<sdsl::bit_vector>();
    delete last_;
    sdsl::bit_vector last(old_last.size() - num_edges_to_remove, false);

    for (edge_index i = 0, new_i = 0; i < edges_to_remove_mask.size(); ++i) {
        if (!edges_to_remove_mask[i]) {
            last[new_i++] = old_last[i];
        } else {
            if (old_last[i] && new_i > 1 && !last[new_i - 1])
                last[new_i - 1] = 1;
        }
    }
    last_ = new bit_vector_stat(std::move(last));

    // update W
    sdsl::int_vector<> old_W = W_->to_vector();
    delete W_;
    sdsl::int_vector<> W(old_W.size() - num_edges_to_remove, 0, bits_per_char_W_);
    std::vector<bool> first_removed(alph_size, false);

    for (edge_index i = 0, new_i = 0; i < edges_to_remove_mask.size(); ++i) {
        TAlphabet c = old_W[i];
        if (edges_to_remove_mask[i]) {
            if (c < alph_size)
                first_removed[c] = true;
        } else {
            assert(c != alph_size);
            if (c > alph_size && first_removed[c % alph_size]) {
                W[new_i++] = c % alph_size;
            } else {
                W[new_i++] = c;
            }
            first_removed[c % alph_size] = false;
        }
    }
    W_ = new wavelet_tree_stat(bits_per_char_W_, std::move(W));

    // update F
    TAlphabet c = 0;
    uint64_t count = 0;
    for (edge_index i = 1; i <= F_.back(); ++i) {
        while (i > F_[c] && c < alph_size) {
            F_[c++] = count;
        }
        if (!edges_to_remove_mask[i])
            count++;
    }
    while (c < alph_size) {
        F_[c++] = count;
    }

    recompute_NF();

    state = State::STAT;

    return num_edges_to_remove;
}

/**
 * Depth first edge traversal.
 * Traverse all edges reachable from the given one.
 */
void BOSS::edge_DFT(edge_index start,
                    Call<edge_index> pre_visit,
                    Call<edge_index> post_visit,
                    std::function<bool(edge_index)> end_branch) const {
    CHECK_INDEX(start);

    // start traversal in the source node of the given edge
    std::vector<edge_index> path { start };
    pre_visit(path.back());

    do {
        // traverse until the last dummy source edge in a path
        while (!end_branch(path.back())) {
            path.push_back(fwd(path.back(), get_W(path.back()) % alph_size));
            pre_visit(path.back());
        }

        // traverse the path backwards to the next branching node
        while (path.size() > 1 && get_last(path.back() - 1)) {
            post_visit(path.back());
            path.pop_back();
        }

        // explore the next edge outgoing from the current branching node
        if (path.size() > 1) {
            post_visit(path.back());
            path.back()--;
            pre_visit(path.back());
        }
    } while (path.size() > 1);

    post_visit(path.back());
}

/**
 * Traverse the entire dummy subtree
 * and find all redundant dummy edges.
 */
template <class Callback>
void traverse_dummy_edges(const BOSS &graph,
                          edge_index subtree_root,
                          size_t check_depth,
                          const Callback &call_edge,
                          sdsl::bit_vector *redundant_mask,
                          sdsl::bit_vector *traversed_mask) {
    assert(!redundant_mask || redundant_mask->size() == graph.get_W().size());
    assert(!traversed_mask || traversed_mask->size() == graph.get_W().size());

    if (!check_depth)
        return;

    constexpr bool async = true;

    std::vector<bool> redundant_path;

    // start traversal in the given node
    graph.edge_DFT(subtree_root,
        [&](edge_index) {
            redundant_path.push_back(true);
        },
        [&](edge_index edge) {
            assert(graph.get_W(edge) < graph.alph_size);

            call_edge(edge, redundant_path.size());

            // TODO: remove this mask and use |call_edge| instead
            if (traversed_mask)
                set_bit(traversed_mask->data(), edge, async, __ATOMIC_RELAXED);

            if (redundant_path.back())
                set_bit(redundant_mask->data(), edge, async);

            redundant_path.pop_back();
        },
        [&](edge_index edge) {
            if (redundant_path.size() == check_depth) {
                if (!redundant_mask
                        || (!fetch_bit(redundant_mask->data(), edge, async)
                                && graph.is_single_incoming(edge, graph.get_W(edge)))) {
                    // the last dummy edge is not redundant and hence the
                    // entire path has to remain in the graph
                    redundant_path.assign(redundant_path.size(), false);
                }
                return true;
            } else {
                assert(graph.is_single_incoming(edge, graph.get_W(edge)));
                return false;
            }
        }
    );
    assert(redundant_path.empty());
}

/**
 * Traverse the entire dummy tree, detect all redundant
 * dummy source edges and return number of these edges.
 */
uint64_t traverse_dummy_edges(const BOSS &graph,
                              sdsl::bit_vector *redundant_mask,
                              sdsl::bit_vector *traversed_mask,
                              size_t num_threads,
                              const std::function<void(edge_index, size_t)> &on_dummy
                                        = [](edge_index, size_t) {}) {
    assert(!redundant_mask || redundant_mask->size() == graph.get_W().size());
    assert(!traversed_mask || traversed_mask->size() == graph.get_W().size());

    // assume the main dummy source node already traversed
    std::atomic<uint64_t> num_dummy_traversed = 1;
    on_dummy(1, 0); // call the main dummy source edge
    if (traversed_mask)
        (*traversed_mask)[1] = true;

    if (graph.get_last(1))
        return 1; // the only dummy k-mer is the main dummy source

    if (graph.get_k() <= 1) {
        // start traversal in the children of the main dummy source
        edge_index root = 2;
        do {
            traverse_dummy_edges(graph, root, graph.get_k(),
                                 [&](edge_index edge, size_t depth) {
                                     num_dummy_traversed++;
                                     on_dummy(edge, depth);
                                 },
                                 redundant_mask, traversed_mask);
            logger->trace("Source dummy edges traversed: {}", num_dummy_traversed);
        } while (!graph.get_last(root++));

        return num_dummy_traversed;
    }

    const size_t tree_split_depth = std::min(size_t(6), graph.get_k() / 2);

    std::vector<edge_index> start_edges;

    // run traversal for subtree at depth |tree_split_depth| in parallel
    edge_index root = 2;
    size_t depth = 0;
    do {
        graph.edge_DFT(root,
            [&](edge_index) { depth++; },
            [&](edge_index) { depth--; },
            [&](edge_index edge) {
                if (depth == tree_split_depth) {
                    start_edges.push_back(edge);
                    return true;
                } else {
                    return false; // go deeper
                }
            }
        );
        assert(!depth);
    } while (!graph.get_last(root++));

    ProgressBar progress_bar(start_edges.size(), "Dummy subtrees traversed",
                             std::cerr, !common::get_verbose());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < start_edges.size(); ++i) {
        traverse_dummy_edges(graph, start_edges[i], 1 + graph.get_k() - tree_split_depth,
            [&](edge_index edge, size_t depth) {
                num_dummy_traversed++;
                on_dummy(edge, tree_split_depth - 1 + depth);
            },
            redundant_mask, traversed_mask
        );
        ++progress_bar;
    }

    logger->trace("All subtrees are traversed");

    // traverse all dummy nodes at depth <= tree_split_depth
    root = 2;
    do {
        traverse_dummy_edges(graph, root, tree_split_depth,
                            [&](edge_index edge, size_t depth) {
                                // We must reach the depth of |tree_split_depth|
                                // to update the masks, but we don't need to
                                // visit those edges, they had already been visited.
                                // So now we visit only their predecessors
                                if (depth < tree_split_depth) {
                                    num_dummy_traversed++;
                                    on_dummy(edge, depth);
                                }
                            },
                            redundant_mask, traversed_mask);
    } while (!graph.get_last(root++));

    return num_dummy_traversed;
}

/**
 * Traverse the entire dummy subgraph (which is a tree)
 * and erase all redundant dummy edges.
 * If passed, mark |source_dummy_edges| with positions
 * of non-redundant dummy source edges.
 * Return value: edges removed from the BOSS graph.
 */
sdsl::bit_vector
BOSS::erase_redundant_dummy_edges(sdsl::bit_vector *source_dummy_edges,
                                  size_t num_threads,
                                  bool verbose) {
    // reset the index of suffix ranges
    index_suffix_ranges(0);

    sdsl::bit_vector redundant_dummy_edges_mask(W_->size(), false);

    if (source_dummy_edges) {
        (*source_dummy_edges) = sdsl::bit_vector();
        (*source_dummy_edges) = sdsl::bit_vector(W_->size(), false);
        (*source_dummy_edges)[1] = true;
    }

    if (get_last(1))
        return redundant_dummy_edges_mask;

    State state = get_state();
    switch_state(State::STAT);

    auto num_dummy_traversed = traverse_dummy_edges(
        *this, &redundant_dummy_edges_mask, source_dummy_edges, num_threads
    );

    if (verbose) {
        std::cout << "Traversal done. Total number of source dummy edges: "
                  << num_dummy_traversed << std::endl;
    }

    auto num_edges_erased = erase_edges(redundant_dummy_edges_mask);
    if (source_dummy_edges)
        utils::erase(source_dummy_edges, redundant_dummy_edges_mask);

    if (verbose) {
        std::cout << "Number of source dummy edges removed: "
                  << num_edges_erased << std::endl;
    }

    switch_state(state);

    return redundant_dummy_edges_mask;
}

uint64_t BOSS::mark_source_dummy_edges(sdsl::bit_vector *mask,
                                       size_t num_threads) const {
    assert(!mask || mask->size() == W_->size());

    return traverse_dummy_edges(*this, NULL, mask, num_threads);
}

uint64_t BOSS::mark_sink_dummy_edges(sdsl::bit_vector *mask) const {
    if (!mask)
        return rank_W(num_edges(), 0) - 1;

    assert(mask->size() == W_->size());

    uint64_t num_dummy_sink_edges = 0;

    // skip the main dummy source
    for (edge_index i = 2; i < W_->size(); ++i) {
        assert(get_W(i) != alph_size);
        if (!get_W(i)) {
            (*mask)[i] = true;
            num_dummy_sink_edges++;
        }
    }

    assert(num_dummy_sink_edges == rank_W(num_edges(), 0) - 1);

    return num_dummy_sink_edges;
}

sdsl::bit_vector BOSS::mark_all_dummy_edges(size_t num_threads) const {
    sdsl::bit_vector edge_mask(num_edges() + 1, 0);

    mark_source_dummy_edges(&edge_mask, num_threads);
    mark_sink_dummy_edges(&edge_mask);

    // exclude 0 as the dummy index that denotes not existing k-mers
    edge_mask[0] = true;

    return edge_mask;
}

sdsl::bit_vector BOSS::prune_and_mark_all_dummy_edges(size_t num_threads) {
    sdsl::bit_vector edge_mask(num_edges() + 1, 0);

    erase_redundant_dummy_edges(&edge_mask, num_threads);
    mark_sink_dummy_edges(&edge_mask);

    // exclude 0 as the dummy index that denotes not existing k-mers
    edge_mask[0] = true;

    return edge_mask;
}

/**
 * Merge BOSS table |other| into the current one. The BOSS |other|
 * is fully traversed and all edges are added to to the current BOSS table.
 * This function is well suited to merge small graphs into large ones.
 */
void BOSS::merge(const BOSS &other, size_t num_threads) {
    std::mutex seq_mutex;
    other.call_sequences([&](const std::string &sequence, auto&&) {
        std::lock_guard<std::mutex> lock(seq_mutex);
        add_sequence(sequence, true);
    }, num_threads);
}

void BOSS::call_start_edges(Call<edge_index> callback) const {
    // start traversal in the main dummy source node and traverse the tree
    // of source dummy k-mers to depth k + 1

    // check if the dummy tree is not empty
    if (is_single_outgoing(1))
        return;

    // run traversal for subtree
    edge_index root = 2;
    size_t depth = 0;
    do {
        edge_DFT(root,
            [&](edge_index) { depth++; },
            [&](edge_index) { depth--; },
            [&](edge_index edge) {
                if (depth < k_)
                    return false;

                if (depth == k_)
                    return !is_single_incoming(edge, get_W(edge));

                callback(edge);
                return true;
            }
        );
        assert(!depth);
    } while (!get_last(root++));
}


// Methods for inferring node degrees with a mask

// If a single outgoing edge is found, write it to |*i| and return true.
// If no outgoing edges are found, set |*i| to 0 and return false.
// If multiple outgoing edges are found, set |*i| to the last and return false.
inline bool masked_pick_single_outgoing(const BOSS &boss,
                                        edge_index *i,
                                        const bitmap *subgraph_mask) {
    assert(i && *i);
    assert(boss.get_last(*i));
    assert(!subgraph_mask || subgraph_mask->size() == boss.num_edges() + 1);

    // in boss, at least one outgoing edge always exists
    if (!subgraph_mask)
        return boss.is_single_outgoing(*i);

    edge_index j = *i;

    bool edge_detected = (*subgraph_mask)[j];

    while (--j > 0 && !boss.get_last(j)) {
        if ((*subgraph_mask)[j]) {
            // stop if there are multiple outgoing edges detected
            if (edge_detected)
                return false;

            // this is the first outgoing edge
            edge_detected = true;
            *i = j;
        }
    }

    // return true of there is exactly one outgoing edge
    if (edge_detected)
        return true;

    // no outgoing
    *i = 0;
    return false;
}

template <class Callback>
inline void masked_call_outgoing(const BOSS &boss,
                                 edge_index i,
                                 const bitmap *subgraph_mask,
                                 const Callback &callback) {
    assert(i);
    assert(boss.get_last(i)
            && "i has to point to the last outgoing edge in unmasked graph");
    assert(!subgraph_mask || subgraph_mask->size() == boss.num_edges() + 1);

    boss.call_outgoing(i, [&](BOSS::edge_index adjacent_edge) {
        if (!subgraph_mask || (*subgraph_mask)[adjacent_edge])
            callback(adjacent_edge);
    });
}

// If a single incoming edge is found, write it to |*i| and return true.
// If no incoming edges are found, set |*i| to 0 and return false.
// If multiple incoming edges are found, set |*i| to the first and return false.
inline bool masked_pick_single_incoming(const BOSS &boss,
                                        edge_index *i, TAlphabet d,
                                        const bitmap *subgraph_mask) {
    assert(i && *i);
    assert(boss.get_W(*i) < boss.alph_size && "must be the first incoming edge");
    assert(d == boss.get_W(*i));
    assert(!subgraph_mask || subgraph_mask->size() == boss.num_edges() + 1);

    // in boss, at least one incoming edge always exists
    if (!subgraph_mask)
        return boss.is_single_incoming(*i, d);

    auto j = *i;

    TAlphabet d_next;
    bool edge_detected = false;
    do {
        if ((*subgraph_mask)[j]) {
            // stop if there are multiple incoming edges detected
            if (edge_detected)
                return false;

            // this is the first incoming edge
            edge_detected = true;
            *i = j;
        }

        if (++j > boss.num_edges())
            break;

        std::tie(j, d_next) = boss.succ_W(j, d, d + boss.alph_size);

    } while (d_next == d + boss.alph_size);

    // return true of there is exactly one incoming edge
    if (edge_detected)
        return true;

    // no incoming edges
    *i = 0;
    return false;
}

using Edge = std::pair<edge_index, std::unique_ptr<TAlphabet[]>>;

// Stores edges during graph traversal, first with decoded k-mers, but
// then switches and only stores edge indices if too many are queued.
class EdgeQueue {
  public:
    EdgeQueue(std::deque<Edge>&& initial_edges = {})
          : decoded_edges_(std::move(initial_edges)) {
        total_decoded_ += decoded_edges_.size();
    }

    explicit EdgeQueue(edge_index edge)
          : indexes_({ edge }) {}

    // don't allow copying
    explicit EdgeQueue(const EdgeQueue &) = delete;
    // allow moving
    EdgeQueue(EdgeQueue&&) = default;

    ~EdgeQueue() {
        total_decoded_ -= decoded_edges_.size();
    }

    template <typename Iterator>
    void emplace_back(edge_index edge, Iterator begin, Iterator end) {
        static_assert(std::is_same_v<std::decay_t<decltype(*begin)>, TAlphabet>);

        if (begin == end) {
            indexes_.push_back(edge);
            return;
        }

        if (++total_decoded_ <= MAX_DECODED_EDGES_IN_QUEUE) {
            decoded_edges_.emplace_back(edge, new TAlphabet[end - begin]);
            std::copy(begin, end, decoded_edges_.back().second.get());
        } else {
            --total_decoded_;
            indexes_.push_back(edge);
        }
    }

    template <typename Iterator>
    void emplace_front(edge_index edge, Iterator begin, Iterator end) {
        static_assert(std::is_same_v<std::decay_t<decltype(*begin)>, TAlphabet>);

        if (begin == end) {
            indexes_.push_front(edge);
            return;
        }

        if (++total_decoded_ <= MAX_DECODED_EDGES_IN_QUEUE) {
            decoded_edges_.emplace_front(edge, new TAlphabet[end - begin]);
            std::copy(begin, end, decoded_edges_.front().second.get());
        } else {
            --total_decoded_;
            indexes_.push_front(edge);
        }
    }

    Edge pop_back() {
        assert(size());
        if (decoded_edges_.size()) {
            Edge edge = std::move(decoded_edges_.back());
            decoded_edges_.pop_back();
            --total_decoded_;
            return edge;
        } else {
            Edge edge(indexes_.back(), nullptr);
            indexes_.pop_back();
            return edge;
        }
    }

    // split off a half of the queue
    // tries to keep in the current queue as many decoded edges as possible
    EdgeQueue split_half() {
        EdgeQueue split_queue;

        // prefer moving edges without k-mers decoded
        size_t old_size = size();
        size_t h = std::min(old_size / 2, indexes_.size());
        split_queue.indexes_.assign(indexes_.end() - h, indexes_.end());
        split_queue.indexes_.shrink_to_fit();
        indexes_.resize(indexes_.size() - h);

        // if moved less than old_size/2 indexes, move decoded k-mers as well
        h = old_size / 2 - h;
        split_queue.decoded_edges_.assign(std::make_move_iterator(decoded_edges_.end() - h),
                                          std::make_move_iterator(decoded_edges_.end()));
        split_queue.decoded_edges_.shrink_to_fit();
        decoded_edges_.resize(decoded_edges_.size() - h);

        return split_queue;
    }

    size_t size() const { return decoded_edges_.size() + indexes_.size(); }

    bool empty() const { return decoded_edges_.empty() && indexes_.empty(); }

  private:
    inline static std::atomic<uint64_t> total_decoded_ = 0;
    std::deque<Edge> decoded_edges_;
    std::deque<edge_index> indexes_;
};

/*
 * Traverse graph from the specified (k+1)-mer/edge and call all paths reachable from it.
 */
void call_paths(const BOSS &boss,
                EdgeQueue&& edges,
                const BOSS::Call<std::vector<edge_index>&&,
                                 std::vector<TAlphabet>&&> &callback,
                bool split_to_unitigs,
                bool kmers_in_single_form,
                bool trim_sentinels,
                ThreadPool &thread_pool,
                sdsl::bit_vector *visited_ptr,
                tsl::hopscotch_set<edge_index> *fetched,
                bool async,
                std::mutex &fetched_mutex,
                ProgressBar &progress_bar,
                const bitmap *subgraph_mask);


void
call_path(const BOSS &boss,
          const BOSS::Call<std::vector<edge_index>&&,
                           std::vector<TAlphabet>&&> &callback,
          std::vector<edge_index> &path,
          std::vector<TAlphabet> &sequence,
          bool split_to_unitigs,
          bool kmers_in_single_form,
          bool trim_sentinels,
          sdsl::bit_vector &visited,
          tsl::hopscotch_set<edge_index> &fetched,
          bool concurrent,
          std::mutex &fetched_mutex,
          ProgressBar &progress_bar,
          const bitmap *subgraph_mask,
          EdgeQueue *edges);

#ifndef NDEBUG
void assert_forks_and_merges_visited(const BOSS &boss,
                                     const sdsl::bit_vector &visited,
                                     const bitmap *subgraph_mask) {
    constexpr bool async = true;
    // make sure that all forks have been covered
    call_zeros(visited, [&](edge_index edge) {
        edge_index t = boss.succ_last(edge);
        bool check = masked_pick_single_outgoing(boss, &t, subgraph_mask);
        assert(t);
        assert(check);

        // make sure the next neighbouring edge has also not been visited
        t = boss.fwd(t, boss.get_W(t) % boss.alph_size);
        check = masked_pick_single_outgoing(boss, &t, subgraph_mask);
        assert(t);
        assert(check);
        assert(!fetch_bit(visited.data(), t, async));
    }, async);

    // make sure that all merges have been covered
    call_zeros(visited, [&](edge_index edge) {
        edge_index t = boss.bwd(edge);
        bool check = masked_pick_single_incoming(boss, &t, boss.get_W(t), subgraph_mask);
        assert(t);
        assert(check);
        assert(!fetch_bit(visited.data(), t, async));
    }, async);
}

void assert_no_leftovers(const BOSS &boss, const sdsl::bit_vector &visited) {
    bool leftover = false;
    call_zeros(visited, [&](edge_index edge) {
        leftover = true;
        TAlphabet d = boss.get_W(edge) % boss.alph_size;
        std::cout << edge << "\t" << boss.get_node_str(edge) << " " << boss.decode(d) << "\t";
        edge = boss.fwd(edge, d);
        d = boss.get_W(edge) % boss.alph_size;
        std::cout << edge << "\t" << boss.get_node_str(edge) << " " << boss.decode(d) << "\n";
    }, true);
    std::cout << std::flush;
    assert(!leftover);
}
#endif

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void BOSS::call_paths(Call<std::vector<edge_index> &&, std::vector<TAlphabet> &&> callback,
                      size_t num_threads,
                      bool split_to_unitigs,
                      bool kmers_in_single_form,
                      const bitmap *subgraph_mask,
                      bool trim_sentinels) const {
    assert(!subgraph_mask || subgraph_mask->size() == W_->size());

    // keep track of the edges that have been reached
    sdsl::bit_vector visited(W_->size(), false);
    if (subgraph_mask) {
        subgraph_mask->add_to(&visited);
        visited.flip();
    }
    visited[0] = true;

    if (trim_sentinels) {
        #pragma omp parallel for num_threads(num_threads) schedule(static, kBlockSize)
        for (uint64_t i = 1; i < W_->size(); ++i) {
            if (get_W(i) == kSentinelCode)
                visited[i] = true;
        }
    }

    tsl::hopscotch_set<edge_index> fetched;
    std::mutex fetched_mutex;

    ProgressBar progress_bar(visited.size() - sdsl::util::cnt_one_bits(visited),
                             "Traverse BOSS",
                             std::cerr, !common::get_verbose());

    ThreadPool thread_pool(num_threads ? num_threads : 1, TASK_POOL_SIZE);
    bool async = true;

    auto enqueue_start = [&](ThreadPool &thread_pool, edge_index start) {
        thread_pool.enqueue([&,start]() {
            ::mtg::graph::boss::call_paths(
                    *this, EdgeQueue(start), callback,
                    split_to_unitigs, kmers_in_single_form,
                    trim_sentinels, thread_pool, &visited, &fetched,
                    async, fetched_mutex, progress_bar, subgraph_mask);
        });
    };

    // start traversal from the source dummy edges first
    //
    //  .____
    //
    // used to indicate if an edge from a source node has been taken
    if (!subgraph_mask) {
        for (edge_index i = succ_last(1); i >= 1; --i) {
            if (!fetch_bit(visited.data(), i, async))
                enqueue_start(thread_pool, i);
        }

    } else {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint64_t begin = 1; begin < visited.size(); begin += kBlockSize) {
            begin = pred_last(begin) + 1;
            // find all edges in a block with no incoming edge
            // and enqueue a traversal task from each one
            uint64_t end = pred_last(std::min(begin + kBlockSize,
                                              visited.size() - 1)) + 1;
            edge_index last_i = 0;
            call_ones(*subgraph_mask, begin, end, [&](edge_index i) {
                i = succ_last(i);
                if (i == last_i)
                    return;

                last_i = i;

                // skip |i| if it has at least one incoming edge
                edge_index j = bwd(i);
                masked_pick_single_incoming(*this, &j, get_W(j), subgraph_mask);
                if (j)
                    return;

                masked_call_outgoing(*this, i, subgraph_mask,
                                     [&](edge_index e) {
                    if (!kmers_in_single_form || !fetch_bit(visited.data(), e, async)) {
                        enqueue_start(thread_pool, e);
                    }
                });
            });
        }
    }

    // then all forks
    //  ____.____
    //       \___
    //
    std::vector<edge_index> edges;
    edge_index last_i = 0;
    call_zeros(visited, [&](edge_index i) {
        if (i <= last_i)
            return;

        last_i = succ_last(i);

        edges.resize(0);
        masked_call_outgoing(*this, last_i, subgraph_mask,
                             [&](edge_index e) { edges.push_back(e); });

        // no outgoing edges or a unique outgoing edge
        if (edges.size() < 2)
            return;

        for (edge_index e : edges) {
            if (!fetch_bit(visited.data(), e, async))
                enqueue_start(thread_pool, e);
        }
    }, async);

    thread_pool.join();

#ifndef NDEBUG
    assert_forks_and_merges_visited(*this, visited, subgraph_mask);
#endif

    // Now we only have to traverse loops that have not been traversed or
    // loops that have partially been traversed (only when extracting
    // primary unitigs/contigs).

    auto process_cycle = [&](edge_index edge) {
        if (fetch_bit(visited.data(), edge, async))
            return;

        edge_index start = edge;
        std::vector<edge_index> path;
        std::vector<TAlphabet> sequence = get_node_seq(edge);
        do {
            TAlphabet w = get_W(edge);
            assert(w != kSentinelCode);
            TAlphabet d = w % alph_size;
            sequence.push_back(d);
            path.push_back(edge);
            edge = fwd(edge, d);
            masked_pick_single_outgoing(*this, &edge, subgraph_mask);
            assert(edge);
        } while (edge != start);

        // Ensures that call_path is called only once for each cycle
        auto rep = std::min_element(path.begin(), path.end());
        if (!fetch_bit(visited.data(), *rep, async)) {
            EdgeQueue queue;
            queue.emplace_back(*rep,
                               sequence.begin() + (rep - path.begin()),
                               sequence.begin() + (rep - path.begin()) + get_k());
            ::mtg::graph::boss::call_paths(
                    *this, std::move(queue), callback,
                    split_to_unitigs, kmers_in_single_form,
                    trim_sentinels, thread_pool, &visited, &fetched,
                    async, fetched_mutex, progress_bar, subgraph_mask);
        }
    };

    if (!async) {
        call_zeros(visited, process_cycle);

    } else {
        std::vector<edge_index> index_buffer;

        call_zeros(visited, [&](edge_index edge) {
            // traverse loops in parallel and only check for unique k-mers at the
            // end of the traversal
            index_buffer.push_back(edge);

            if (index_buffer.size() == TRAVERSAL_START_BATCH_SIZE) {
                thread_pool.enqueue([&,index_buffer]() {
                    std::for_each(index_buffer.begin(), index_buffer.end(),
                                  process_cycle);
                });

                index_buffer.clear();
            }

        }, async);

        thread_pool.enqueue([&,index_buffer]() {
            std::for_each(index_buffer.begin(), index_buffer.end(),
                          process_cycle);
        });

        thread_pool.join();
    }

#ifndef NDEBUG
    assert_no_leftovers(*this, visited);
#endif
}

void call_paths(const BOSS &boss,
                EdgeQueue&& edges,
                const BOSS::Call<std::vector<edge_index>&&,
                                 std::vector<TAlphabet>&&> &callback,
                bool split_to_unitigs,
                bool kmers_in_single_form,
                bool trim_sentinels,
                ThreadPool &thread_pool,
                sdsl::bit_vector *visited_ptr,
                tsl::hopscotch_set<edge_index> *fetched_ptr,
                bool async,
                std::mutex &fetched_mutex,
                ProgressBar &progress_bar,
                const bitmap *subgraph_mask) {
    assert(visited_ptr && fetched_ptr);

    auto &visited = *visited_ptr;
    // store all branch nodes on the way
    std::vector<edge_index> out_edges;

    // keep traversing until we have worked off all branches from the queue
    while (!edges.empty()) {
        auto [edge, kmer] = edges.pop_back();

        if (fetch_bit(visited.data(), edge, async))
            continue;

        std::vector<edge_index> path;
        path.reserve(100);

        std::vector<TAlphabet> sequence;
        if (kmer) {
            sequence.reserve(100 + boss.get_k());
            sequence.assign(kmer.get(), kmer.get() + boss.get_k());
            kmer.reset();
        } else {
            sequence = boss.get_node_seq(edge);
            sequence.reserve(100 + boss.get_k());
        }

        // traverse simple path until we reach its tail or
        // the first edge that has already been visited
        while (!fetch_and_set_bit(visited.data(), edge, async)) {
            assert(edge > 0);
            assert(!subgraph_mask || (*subgraph_mask)[edge]);
            ++progress_bar;

            // visit the edge
            TAlphabet w = boss.get_W(edge);
            TAlphabet d = w % boss.alph_size;
            sequence.push_back(d);
            path.push_back(edge);

            // stop the traversal if the next node is a dummy sink
            if (!d)
                break;

            bool stop_even_if_single_outgoing;

            if (!split_to_unitigs) {
                // we always continue traversal of a contig if there is
                // only a single outgoing edge
                stop_even_if_single_outgoing = false;

            } else if (!subgraph_mask && w != d) {
                // For unitigs, we must terminate the traversal if there
                // are more than one edges incoming to the target node.
                // If the entire graph is selected and the edge is marked
                // with '-', we know for sure this is not the only edge
                // incoming into its target node.
                stop_even_if_single_outgoing = true;
            } else {
                // Otherwise, we must check all edges by iteration.
                // shift to the first incoming edge
                if (w != d)
                    edge = boss.pred_W(edge, d);
                // check if there are multiple incoming edges
                stop_even_if_single_outgoing
                    = !masked_pick_single_incoming(boss, &edge, d, subgraph_mask);
                assert(edge);
            }

            // make one traversal step
            edge = boss.fwd(edge, d);

            out_edges.resize(0);
            masked_call_outgoing(boss, edge, subgraph_mask,
                                 [&](edge_index e) { out_edges.push_back(e); });
            // stop the traversal if there are no edges outgoing from the target
            if (out_edges.empty())
                break;

            edge = out_edges.front();

            // continue the non-branching traversal further if
            //      1. there is only one edge outgoing from the target
            // and at least one of the following conditions is met:
            //      2. there is only one edge incoming to the target node
            //      3. we call contigs (the unitigs may be concatenated)
            if (out_edges.size() == 1 && !stop_even_if_single_outgoing)
                continue;

            edge_index next_edge = 0;

            // loop over the outgoing edges
            for (edge_index edge : out_edges) {
                assert((!subgraph_mask || (*subgraph_mask)[edge]
                    || fetch_bit(visited.data(), edge, async))
                        && "k-mers not from subgraph are marked as visited");

                if (!fetch_bit(visited.data(), edge, async)) {
                    if (!next_edge && !split_to_unitigs) {
                        // save the edge for visiting if we extract contigs
                        next_edge = edge;
                    } else {
                        // visit other edges
                        edges.emplace_back(edge, sequence.end() - boss.get_k(), sequence.end());
                    }
                }
            }

            // stop traversing this sequence if the next edge was not selected
            if (!next_edge)
                break;

            // pick the last outgoing but not yet visited
            // edge and continue traversing the graph
            edge = next_edge;

            if (edges.size() >= MAX_EDGE_QUEUE_SIZE) {
                thread_pool.force_enqueue(
                    [=,&boss,&thread_pool,&fetched_mutex,&progress_bar](EdgeQueue &edges) {
                        ::mtg::graph::boss::call_paths(
                                boss, std::move(edges), callback,
                                split_to_unitigs, kmers_in_single_form,
                                trim_sentinels, thread_pool, visited_ptr, fetched_ptr,
                                async, fetched_mutex, progress_bar, subgraph_mask);
                    },
                    edges.split_half()
                );
            }
        }

        if (path.empty())
            continue;

        call_path(boss, callback, path, sequence, split_to_unitigs,
                  kmers_in_single_form, trim_sentinels, visited, *fetched_ptr,
                  async, fetched_mutex, progress_bar, subgraph_mask, &edges);
    }
}

/**
 * Updates #terminal and #near_terminal based on the given path.
 * One terminal node is set every max_length nodes in the path, and all nodes before it
 * are marked in #near_terminal.
 * The last node in the path is marked as terminal if:
 *  1. The path length is an exact multiple of max_length, OR
 *  2. The last node in the path merges into a node that is neither terminal nor near
 *     terminal
 */
void update_terminal_bits(size_t max_length,
                          edge_index next_edge,
                          std::vector<edge_index> &&path,
                          sdsl::bit_vector *terminal,
                          sdsl::bit_vector *near_terminal) {
    assert(next_edge);

    if (path.empty())
        return;

    size_t i = 0;
    constexpr bool async = true;

    // set anchors
    // .........V..........V....*
    // ||||||||||
    // max_length
    for ( ; i + max_length <= path.size(); i += max_length) {
        for (uint64_t j = i; j + 1 < i + max_length; ++j) {
            assert(!fetch_bit(terminal->data(), path[j], async));
            assert(!fetch_bit(near_terminal->data(), path[j], async));
            set_bit(near_terminal->data(), path[j], async);
        }
        assert(!fetch_bit(near_terminal->data(), path[i + max_length - 1], async));
        set_bit(terminal->data(), path[i + max_length - 1], async);
    }

    if (i == path.size()) // last node is terminal
        return;

    // current position |i|
    // .........V..........V....*
    //                      ^   ^
    //                      i next_edge

    // If the next node is close to an existing anchor and, therefore, every
    // node in this last segment is close to it too (at most 2 * max_length),
    // there is no need to set another anchor at the end of the path.
    if (fetch_bit(near_terminal->data(), next_edge, async)) {
        assert(!fetch_bit(terminal->data(), next_edge, async));
        return;
    }

    if (!fetch_bit(terminal->data(), next_edge, async)) {
        set_bit(terminal->data(), path.back(), async);
        path.pop_back();
    }

    // The last node of |path| is an anchor, so all its predecessors must be
    // marked as near-anchor.
    while (i < path.size()) {
        set_bit(near_terminal->data(), path[i++], async);
    }
}

/**
 * Traverses the path that can be visited starting from #edge.
 * A path ends when there are either no outgoing edges from the current node or
 * if the row-diff successor in a fork has already been visited.
 */
void traverse_rd_path(const BOSS &boss,
                      const bit_vector &rd_succ,
                      edge_index edge,
                      size_t max_length,
                      sdsl::bit_vector *visited,
                      sdsl::bit_vector *terminal,
                      sdsl::bit_vector *near_terminal,
                      ProgressBar &progress_bar) {
    assert(visited && terminal && near_terminal);

    // make sure it's not a dummy k-mer
    assert(boss.get_node_seq(edge).front() != boss.kSentinelCode);

    constexpr bool async = true;

    if (fetch_bit(visited->data(), edge, async))
        return;

    std::vector<edge_index> path;
    path.reserve(100);

    // traverse unvisited simple path, always pick a successor according
    // to the routing in row-diff
    while (!fetch_and_set_bit(visited->data(), edge, async)) {
        assert(edge > 0);
        ++progress_bar;

        path.push_back(edge);

        edge = boss.row_diff_successor(edge, rd_succ);
    }

    // mark terminal and near terminal nodes
    update_terminal_bits(max_length, edge, std::move(path), terminal, near_terminal);
}

// Returns new edges visited while fetching the path (only returns
// a non-empty set for primary mode |kmers_in_single_form| = true).
// Since fwd will be called on all edges in the returned vector, the corresponding
// node sequences have been precomputed in these Edges
// e.g.,
// edge.first: ATGGGT G -> edge.second = {T,G,G,G,T,G}
void
call_path(const BOSS &boss,
          const BOSS::Call<std::vector<edge_index>&&,
                           std::vector<TAlphabet>&&> &callback,
          std::vector<edge_index> &path,
          std::vector<TAlphabet> &sequence,
          bool split_to_unitigs,
          bool kmers_in_single_form,
          bool trim_sentinels,
          sdsl::bit_vector &visited,
          tsl::hopscotch_set<edge_index> &fetched,
          bool concurrent,
          std::mutex &fetched_mutex,
          ProgressBar &progress_bar,
          const bitmap *subgraph_mask,
          EdgeQueue *edge_queue) {
#ifndef NDEBUG
    for (edge_index e : path) {
        assert(e);
        assert(fetch_bit(visited.data(), e, concurrent));
    }
#endif

    if (!trim_sentinels && !kmers_in_single_form) {
        callback(std::move(path), std::move(sequence));
        return;
    }

    // trim trailing sentinels '$'
    if (sequence.back() == boss.kSentinelCode) {
        sequence.pop_back();
        path.pop_back();
    }

    auto first_valid_it
        = std::find_if(sequence.begin(), sequence.end(),
                       [&boss](auto c) { return c != boss.kSentinelCode; });

    if (first_valid_it + boss.get_k() >= sequence.end())
        return;

    sequence.erase(sequence.begin(), first_valid_it);
    path.erase(path.begin(),
               path.begin() + (first_valid_it - sequence.begin()));

    if (!kmers_in_single_form) {
        callback(std::move(path), std::move(sequence));
        return;
    }

    // get dual path (mapping of the reverse complement sequence)
    auto rev_comp_seq = sequence;
    kmer::KmerExtractorBOSS::reverse_complement(&rev_comp_seq);

    auto dual_path = boss.map_to_edges(rev_comp_seq);
    // restrict to the subgraph
    for (edge_index &e : dual_path) {
        if (subgraph_mask && !(*subgraph_mask)[e]) {
            e = 0;
        }
    }

    // stores the dual nodes that had been visited by other threads, and hence
    // the following contention must be resolved:
    // this thread: primary vs dual,
    // other thread: dual vs primary.
    std::vector<edge_index> dual_visited;
    dual_visited.reserve(path.size());

    // first, we mark all reverse-complement (dual) k-mers as visited
    for (size_t i = 0; i < dual_path.size(); ++i) {
        if (!dual_path[i])
            continue;

        if (!fetch_and_set_bit(visited.data(), dual_path[i], concurrent)) {
            ++progress_bar;

            // schedule traversal branched off from each terminal dual k-mer
            if (i + 1 == dual_path.size() || !dual_path[i + 1]) {
                edge_index next_edge = boss.fwd(dual_path[i],
                                                boss.get_W(dual_path[i]) % boss.alph_size);

                // schedule only if it has a single outgoing k-mer,
                // otherwise it's a fork which will be covered in the forward pass.
                if (masked_pick_single_outgoing(boss, &next_edge, subgraph_mask)) {
                    if (!fetch_bit(visited.data(), next_edge, concurrent)) {
                        edge_queue->emplace_front(next_edge,
                                                 rev_comp_seq.begin() + i + 1,
                                                 rev_comp_seq.begin() + i + 1 + boss.get_k());
                    }
                }
            }
        } else {
            // The dual node had already been visited, so we insert its index
            // to the buffer. The index is inverted because dual_path will be
            // reversed.
            dual_visited.push_back(dual_path.size() - 1 - i);

            // For the unitig mode, wait until the dual unitig is fully traversed
            // by the same thread that visited its first node.
            if (i == 0 && split_to_unitigs
                       && !std::count(dual_path.begin(), dual_path.end(), 0)) {
                // The first node had already been visited, hence, the remaining
                // part of this non-branching path will be reached as well.
                while (!fetch_bit(visited.data(), dual_path.back(), concurrent)) {}
            }
        }
    }

    if (!dual_visited.size()) {
        callback(std::move(path), std::move(sequence));
        return;
    }

    std::reverse(dual_path.begin(), dual_path.end());

    // find all the fetched points where the path must be cut
    std::vector<size_t> breakpoints;
    breakpoints.reserve(dual_visited.size());

    // then lock all threads and resolve all competing nodes from dual_visited
    {
        std::unique_lock<std::mutex> lock(fetched_mutex);

        for (size_t i : dual_visited) {
            if (!fetched.count(dual_path[i])) {
                assert(!fetched.count(path[i]));
                // the dual node had not been fetched by the other thread, so
                // we'll fetch path[i] in this thread and mark it as visited,
                // so the other thread will find it in the hash set and erase
                // it from there
                fetched.insert(path[i]);

            } else {
                // the dual node had been fetched first by the other thread, so
                // we'll skip path[i] in this thread
                breakpoints.push_back(i);
                // erase the node from the hash set
                fetched.erase(dual_path[i]);
            }
        }
    }

    // sort breakpoints (the initial order was derived from the dual path, hence reversed)
    std::reverse(breakpoints.begin(), breakpoints.end());
    // include the last segment
    breakpoints.push_back(path.size());

    // fetch the segments cut off from the path
    size_t begin = 0;
    for (size_t i : breakpoints) {
        // The k-mer or its reverse-complement k-mer had been fetched
        // -> Skip this k-mer and call the traversed path segment.
        if (begin < i) {
            callback({ path.begin() + begin, path.begin() + i },
                     { sequence.begin() + begin, sequence.begin() + i + boss.get_k() });
        }

        begin = i + 1;
    }
}

void BOSS::call_sequences(Call<std::string&&, std::vector<edge_index>&&> callback,
                          size_t num_threads,
                          bool kmers_in_single_form,
                          const bitmap *subgraph_mask) const {
    call_paths([&](std::vector<edge_index>&& edges, std::vector<TAlphabet>&& path) {
        assert(path.size() >= k_ + 1);
        assert(edges.size() == path.size() - k_);
        assert(!std::count(path.begin(), path.end(), kSentinelCode));

        std::string sequence(path.size(), '\0');
        std::transform(path.begin(), path.end(), sequence.begin(),
                       [&](TAlphabet c) { return BOSS::decode(c); });

        callback(std::move(sequence), std::move(edges));

    }, num_threads, false, kmers_in_single_form, subgraph_mask, true);
}

// Reach all k-mers that merge into anchor |edge| by following their diff paths.
template <typename T>
void traverse_rd_path_backward(const BOSS &boss,
                               const bit_vector &rd_succ,
                               edge_index edge,
                               T max_length,
                               sdsl::bit_vector *visited,
                               sdsl::bit_vector *terminal,
                               sdsl::bit_vector *dummy,
                               ProgressBar &progress_bar) {
    constexpr bool async = true;

    assert(max_length);
    assert(fetch_bit(terminal->data(), edge, async));

    std::vector<std::pair<edge_index, T /* next anchor dist */>> queue;
    queue.emplace_back(edge, 0);

    while (queue.size()) {
        T dist_to_anchor;
        std::tie(edge, dist_to_anchor) = queue.back();
        queue.pop_back();

        assert(boss.get_W(edge) && !fetch_bit(dummy->data(), edge, async));

        // mark as visited
        // also check if the node had already been visited (in case of loop)
        if (fetch_and_set_bit(visited->data(), edge, async))
            continue;

        ++progress_bar;

        if (dist_to_anchor == max_length) {
            // make this node an anchor
            set_bit(terminal->data(), edge, async);
            dist_to_anchor = 0;
        } else if (fetch_bit(terminal->data(), edge, async)) {
            dist_to_anchor = 0;
        }

        // stop if the edge is not the row-diff successor of its source node
        // AAAX - AAX$
        // ^^^^
        // AAAY - ****
        if (!rd_succ[edge])
            continue;

        // |edge| is the row-diff successor. Thus, it is part of a diff
        // path, and all edges incoming to it will compute diff wrt to it.
        // So, we propagate the diff path backward.
        boss.call_incoming_to_target(boss.bwd(edge), boss.get_node_last_value(edge),
            [&](edge_index pred) {
                if (!fetch_bit(dummy->data(), pred, async))
                    queue.emplace_back(pred, dist_to_anchor + 1);
            }
        );
    }
}

void BOSS::row_diff_traverse(size_t num_threads,
                             size_t max_length,
                             const bit_vector &rd_succ,
                             sdsl::bit_vector *terminal) const {
    // TODO: can we do it without using this extra |dummy| bitmap?
    sdsl::bit_vector dummy(W_->size(), false);
    dummy[0] = true;
    // keep track of the edges that have been reached
    sdsl::bit_vector visited = dummy;

    constexpr bool async = true;

    // traverse all dummy source k-mers but leave last source dummy edges
    // as not visited
    traverse_dummy_edges(*this, NULL, NULL, num_threads,
        [&](edge_index edge, size_t depth) {
            assert(depth <= get_k());
            set_bit(dummy.data(), edge, async);
            if (depth < get_k())
                set_bit(visited.data(), edge, async);
        }
    );

    assert(terminal->size() == W_->size());

    // TODO: can we use |near_terminal| to mark visited?
    sdsl::bit_vector near_terminal(W_->size(), false);

    ProgressBar progress_bar(visited.size() - sdsl::util::cnt_one_bits(visited),
                             "Traverse graph", std::cerr, !common::get_verbose());

    // anchor all sinks (incoming to X...X$)
    //  ____.
    #pragma omp parallel for num_threads(num_threads)
    for (edge_index i = 2; i < W_->size(); ++i) {
        if (!get_W(i)) {
            // mark the dummy sink
            set_bit(dummy.data(), i, async);

            assert(!fetch_bit(visited.data(), i, async));
            // mark as visited
            set_bit(visited.data(), i, async);
            ++progress_bar;

            call_incoming_to_target(bwd(i), get_node_last_value(i),
                [&](edge_index pred) {
                    assert(!fetch_bit(dummy.data(), pred, async));
                    // anchor the sink
                    set_bit(terminal->data(), pred, async);
                }
            );
        }
    }

    ThreadPool thread_pool(std::max(num_threads, 1UL), TASK_POOL_SIZE);

    std::function<void(edge_index)> traverse_path;
    std::vector<uint64_t> to_visit;
    auto flush_batch = [&]() {
        thread_pool.enqueue([to_visit,&traverse_path]() {
            std::for_each(to_visit.begin(), to_visit.end(), traverse_path);
        });
        to_visit.resize(0);
    };
    auto enqueue_start = [&](edge_index start) {
        to_visit.push_back(start);
        if (to_visit.size() == TRAVERSAL_START_BATCH_SIZE)
            flush_batch();
    };

    // backward traversal
    traverse_path = [&](edge_index anchor) {
        traverse_rd_path_backward(*this, rd_succ, anchor, max_length, &visited,
                                  terminal, &dummy, progress_bar);
    };

    // run backward traversal from every anchor
    call_ones(*terminal, enqueue_start, async);

    flush_batch();
    thread_pool.join();

    // forward traversal
    traverse_path = [&](edge_index start) {
        traverse_rd_path(*this, rd_succ, start, max_length, &visited,
                         terminal, &near_terminal, progress_bar);
    };
    // start traversal from the dummy source edges first ($X...X)
    // they are marked as |dummy| AND NOT |visited|
    //  .____
    call_ones(dummy, [&](edge_index i) {
        if (!fetch_and_set_bit(visited.data(), i, async)) {
            ++progress_bar;
            edge_index real_source = fwd(i, get_W(i) % alph_size);
            masked_call_outgoing(*this, real_source, NULL, [&](edge_index next) {
                if (!fetch_bit(visited.data(), next, async))
                    enqueue_start(next);
            });
        }
    }, async);

    // start traversal from the last anchors
    //  .....V---->
    call_ones(*terminal, [&](edge_index i) {
        assert(fetch_bit(visited.data(), i, async));
        if (TAlphabet d = get_W(i) % alph_size) {
            masked_call_outgoing(*this, fwd(i, d), NULL, [&](edge_index next) {
                if (!fetch_bit(visited.data(), next, async))
                    enqueue_start(next);
            });
        }
    }, async);

    // then all forks
    //  ____.____
    //       \___
    uint64_t last_processed = 0;
    call_zeros(visited, [&](edge_index i) {
        if (i <= last_processed)
            return; // this fork was already processed

        last_processed = succ_last(i);
        if (last_processed < 2 || get_last(last_processed - 1)) {
            return; // single outgoing edge, so not a fork
        }
        for (; i <= last_processed; ++i) {
            if (!fetch_bit(visited.data(), i, async))
                enqueue_start(i);
        }
    }, async);

    flush_batch();
    thread_pool.join();

#ifndef NDEBUG
    assert_forks_and_merges_visited(*this, visited, nullptr);
#endif

    // Now we only have to traverse simple cycles that have no forks
    traverse_path = [&](edge_index edge) {
        if (fetch_bit(visited.data(), edge, async))
            return;

        edge_index start = edge;
        std::vector<edge_index> path;
        std::vector<TAlphabet> sequence = get_node_seq(edge);
        do {
            TAlphabet w = get_W(edge);
            assert(w != kSentinelCode);
            TAlphabet d = w % alph_size;
            sequence.push_back(d);
            path.push_back(edge);
            assert(is_single_outgoing(edge));
            edge = fwd(edge, d);
        } while (edge != start);

        // check the cycle's representative node to see if the cycle has already been
        // visited
        edge_index rep = *std::min_element(path.begin(), path.end());
        if (fetch_and_set_bit(visited.data(), rep, async))
            return;

        for (edge_index idx : path) {
            std::ignore = idx;
            assert(idx == rep || !fetch_bit(visited.data(), idx, async));
            set_bit(visited.data(), idx, async);
        }

        progress_bar += path.size();

        for (uint64_t i = 0; i < path.size(); i += max_length) {
            set_bit(terminal->data(), path[i], async);
        }
    };

    // traverse cycles in parallel
    call_zeros(visited, enqueue_start, async);

    flush_batch();
    thread_pool.join();

#ifndef NDEBUG
    assert_no_leftovers(*this, visited);
#endif
}

void BOSS::call_unitigs(Call<std::string&&, std::vector<edge_index>&&> callback,
                        size_t num_threads,
                        size_t min_tip_size,
                        bool kmers_in_single_form,
                        const bitmap *subgraph_mask) const {
    call_paths([&](std::vector<edge_index>&& edges, std::vector<TAlphabet>&& path) {
        assert(path.size() >= k_ + 1);
        assert(edges.size() == path.size() - k_);
        assert(!std::count(path.begin(), path.end(), kSentinelCode));

        std::string sequence(path.size(), '\0');
        std::transform(path.begin(), path.end(), sequence.begin(),
                       [&](TAlphabet c) { return BOSS::decode(c); });

        auto first_edge = edges.front();
        auto last_edge = edges.back();

        // always call long unitigs
        if (sequence.size() >= k_ + min_tip_size) {
            callback(std::move(sequence), std::move(edges));
            return;
        }

        // this is a short unitig that must be skipped if it is a tip

        /**
         * Paths are sequences of edges!
         *
         * dead end:
         *          _._._._$
         *
         * both 1 and 2 are dead ends:
         *          1 ..._._._$
         *          2 ..._._./
         *
         * both 1 and 2 are dead ends:
         *          1 ..._._._._$
         *          2 ..._._./
         *
         * neither 1 nor 2 is a dead end:
         *          1 ..._._._._._$
         *          2 ..._._./
         */

        edge_index last_fwd = 0;

        // if the last node has multiple outgoing edges,
        // it is clearly neither a sink tip nor a source tip.
        if (path.back() != kSentinelCode
                && !masked_pick_single_outgoing(*this,
                                                &(last_fwd = fwd(last_edge, path.back())),
                                                subgraph_mask)
                && last_fwd) {
            callback(std::move(sequence), std::move(edges));
            return;
        }

        edge_index first_bwd = 0;

        // if the first node has multiple incoming edges,
        // it is clearly neither a source tip nor a sink tip.
        // TODO: This doesn't work properly if graph has redundant dummy edges.
        //       Make sure there are no redundant dummy edges when this
        //       function is called.
        if (path.front() != kSentinelCode
                && !masked_pick_single_incoming(*this,
                                                &(first_bwd = bwd(first_edge)),
                                                get_node_last_value(first_edge),
                                                subgraph_mask)
                && first_bwd) {
            callback(std::move(sequence), std::move(edges));
            return;
        }

        // this unitig has only one incoming and one
        // outgoing edges (which may be dummy edges)

        // skip all sink dead ends, as they are also sink
        // tips (because there is only one edge incoming to the first node)
        if (path.back() == kSentinelCode
                || !last_fwd
                || !get_W(last_fwd))
            return;

        // skip all source dead ends, as they are also source
        // tips (because there is only one edge outgoing from the last node)
        if (path.front() == kSentinelCode
                || !first_bwd
                || !get_minus_k_value(first_bwd, k_ - 1).first)
            return;

        // this is not a tip
        callback(std::move(sequence), std::move(edges));

    }, num_threads, true, kmers_in_single_form, subgraph_mask, true);
}

/**
 * Traverse the boss graph and call all its edges
 * except for the dummy source nodes and the dummy sink nodes
 */
void BOSS::call_kmers(Call<edge_index, const std::string&> callback) const {
    sdsl::bit_vector visited(W_->size(), false);

    // store all branch nodes on the way
    std::queue<std::pair<edge_index, std::string>> branchnodes;

    // start from the second edge (skip dummy main source)
    for (edge_index i = 2; i < W_->size(); ++i) {
        if (visited[i] || !get_last(i))
            continue;

        //TODO: traverse backwards

        // TODO: bound size of queue to reduce memory usage
        branchnodes.push({ i, get_node_str(i) + '\0' });

        // keep traversing until we have worked off all branches from the queue
        while (!branchnodes.empty()) {
            auto [edge, kmer] = std::move(branchnodes.front());
            branchnodes.pop();

            // traverse forwards until we reach a sink or
            // the first edge that has been already visited
            while (!visited[edge]) {
                assert(edge > 0);
                assert(get_last(edge));

                visited[edge] = true;

                TAlphabet d = get_W(edge) % alph_size;

                // stop traversing if it's a sink
                if (!d)
                    break;

                // traverse if there is only one outgoing edge
                if (is_single_outgoing(edge)) {
                    auto next_edge = fwd(edge, d);

                    kmer.back() = decode(get_node_last_value(next_edge));
                    if (kmer.front() != BOSS::kSentinel)
                        callback(edge, kmer);

                    edge = next_edge;
                    std::copy(kmer.begin() + 1, kmer.end(), kmer.begin());

                } else {
                    // loop over outgoing edges
                    do {
                        assert(get_W(edge));

                        auto next_edge = fwd(edge, d);

                        kmer.back() = decode(get_node_last_value(next_edge));
                        if (kmer.front() != BOSS::kSentinel)
                            callback(edge, kmer);

                        if (!visited[next_edge])
                            branchnodes.push({ next_edge, kmer.substr(1) + '\0' });

                    } while (--edge > 1 && !get_last(edge) && (d = get_W(edge) % alph_size));

                    break;
                }
            }
        }
    }
}

void BOSS::index_suffix_ranges(size_t suffix_length) {
    assert(suffix_length <= k_);

    indexed_suffix_length_ = suffix_length;
    indexed_suffix_ranges_.clear();

    if (indexed_suffix_length_ == 0u)
        return;

    if (indexed_suffix_length_ * log2(alph_size - 1) >= 64)
        throw std::runtime_error("ERROR: Trying to index too long suffixes");

    std::vector<std::tuple<uint64_t, edge_index, edge_index>> suffix_ranges;

    // first, take empty suffix and the entire range of nodes in the BOSS table
    uint64_t num_suffixes = 1;
    suffix_ranges.emplace_back(0, 1, W_->size() - 1);

    // grow the suffix length up to |suffix_length - 1|
    for (size_t len = 1; len < indexed_suffix_length_; ++len) {
        std::vector<std::tuple<uint64_t, edge_index, edge_index>> narrowed;
        narrowed.reserve(suffix_ranges.size() * (alph_size - 1));

        for (const auto &[idx, rl, ru] : suffix_ranges) {
            // prepend the suffix with one of the |alph_size - 1| possible characters
            for (TAlphabet c = 1; c < alph_size; ++c) {
                edge_index rl_next = rl;
                edge_index ru_next = ru;
                if (!tighten_range(&rl_next, &ru_next, c))
                    continue;

                assert(idx < num_suffixes);
                narrowed.emplace_back(num_suffixes * (c - 1) + idx,
                                      rl_next, ru_next);
            }
        }
        suffix_ranges.swap(narrowed);
        num_suffixes *= (alph_size - 1);
    }

    // grow the suffix length the last time and build the final index
    indexed_suffix_ranges_.assign(num_suffixes * (alph_size - 1),
                                  std::pair<edge_index, edge_index>(W_->size(), 0));

    for (const auto &[idx, rl, ru] : suffix_ranges) {
        // prepend the suffix with one of the |alph_size - 1| possible characters
        for (TAlphabet c = 1; c < alph_size; ++c) {
            // tighten the range
            edge_index rl_next = rl;
            edge_index ru_next = ru;
            if (!tighten_range(&rl_next, &ru_next, c))
                continue;

            assert(idx < num_suffixes);
            indexed_suffix_ranges_[num_suffixes * (c - 1) + idx]
                = std::make_pair(rl_next, ru_next);
        }
    }

    // align the upper bounds to enable the binary search on them
    for (size_t i = 1; i < indexed_suffix_ranges_.size(); ++i) {
        if (!indexed_suffix_ranges_[i].second) {
            // shift the upper bounds of the empty ranges but still keep them empty
            indexed_suffix_ranges_[i].second = indexed_suffix_ranges_[i - 1].second;
            indexed_suffix_ranges_[i].first = indexed_suffix_ranges_[i - 1].second + 1;
        }
    }
    assert(std::is_sorted(indexed_suffix_ranges_.begin(),
                          indexed_suffix_ranges_.end(),
                          utils::LessSecond()));
}

bool BOSS::is_valid() const {
    assert((*W_)[0] == 0);
    assert(W_->size() >= 2);
    assert(get_node_str(1) == std::string(k_, kSentinel) && "First kmer must be dummy");
    assert(get_W(1) == kSentinelCode && "First kmer must be dummy");

    for (edge_index i = 1; i < W_->size(); i++) {
        if (get_node_last_value(i) >= alph_size
                || get_W(i) == alph_size
                || get_W(i) >= 2 * alph_size)
            return false;

        auto index_pred = bwd(i);
        if (index_pred < 1
                || index_pred >= W_->size()
                || get_node_last_value(index_pred) >= alph_size
                || get_W(index_pred) == alph_size
                || get_W(index_pred) >= 2 * alph_size)
            return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream &os, const BOSS &graph) {
    graph.print(os);
    return os;
}

} // namespace boss
} // namespace graph
} // namespace mtg
