#include "boss.hpp"

#include <cassert>
#include <vector>
#include <stack>
#include <algorithm>
#include <string>
#include <cstdio>

#include <progress_bar.hpp>
#include <boost/multiprecision/integer.hpp>
#include <libmaus2/util/NumberSerialisation.hpp>

#include "boss_construct.hpp"
#include "serialization.hpp"
#include "reverse_complement.hpp"
#include "utils.hpp"
#include "threading.hpp"

using utils::remove_suffix;
using TAlphabet = BOSS::TAlphabet;


#define CHECK_INDEX(idx) \
    assert(idx < W_->size()); \
    assert(idx > 0)
#define CHECK_NODE(idx) \
    assert(idx <= num_nodes()); \
    assert(idx > 0)

typedef BOSS::node_index node_index;
typedef BOSS::edge_index edge_index;

// TODO: run benchmarks and optimize these parameters
const size_t MAX_ITER_WAVELET_TREE_STAT = 1000;
const size_t MAX_ITER_WAVELET_TREE_DYN = 0;
const size_t MAX_ITER_WAVELET_TREE_SMALL = 10;


BOSS::BOSS(size_t k)
      : alph_size(kmer_extractor_.alphabet.size()),
        alphabet(kmer_extractor_.alphabet),
        bits_per_char_W_(boost::multiprecision::msb(alph_size - 1) + 2),
        k_(k),
        last_(new bit_vector_dyn()),
        F_(alph_size, 0),
        W_(new wavelet_tree_dyn(bits_per_char_W_)) {

    assert(bits_per_char_W_ <= sizeof(TAlphabet) * 8
            && "Choose type for TAlphabet properly");

    last_->insert_bit(0, false);
    W_->insert(0, 0);

    // add the dummy source node
    last_->insert_bit(1, true);
    W_->insert(0, 0);
    for (size_t j = 1; j < alph_size; j++) {
        F_[j] = 1;
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
    for (uint64_t i = 0; i < W_->size(); ++i) {
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
    for (uint64_t i = 0; i < W_->size(); ++i) {
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
    for (uint64_t i = 0; i < F_.size(); ++i) {
        if (get_F(i) != other.get_F(i)) {
            if (verbose)
                std::cout << "F differs at position " << i
                          << "\n1: F[" << i << "] = " << get_F(i)
                          << "\n2: F[" << i << "] = " << other.get_F(i)
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
    uint64_t i = 1;
    uint64_t j = 1;

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
    const auto out_filename = remove_suffix(filename, kExtension) + kExtension;

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
    libmaus2::util::NumberSerialisation::serialiseNumberVector(outstream, F_);
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
    auto file = remove_suffix(filename, kExtension) + kExtension;

    std::ifstream instream(file, std::ios::binary);

    return load(instream);
}

bool BOSS::load(std::ifstream &instream) {
    // if not specified in the file, the default for loading is dynamic
    state = Config::DYN;

    try {
        // load F, k, and state
        F_ = libmaus2::util::NumberSerialisation::deserialiseNumberVector<uint64_t>(instream);
        k_ = load_number(instream);
        state = static_cast<Config::StateType>(load_number(instream));

        if (F_.size() != alph_size) {
            std::cerr << "ERROR: failed to load F vector, incompatible size" << std::endl;
            return false;
        }

        // load W and last arrays
        delete W_;
        delete last_;
        switch (state) {
            case Config::DYN:
                W_ = new wavelet_tree_dyn(bits_per_char_W_);
                last_ = new bit_vector_dyn();
                break;
            case Config::STAT:
                W_ = new wavelet_tree_stat(bits_per_char_W_);
                last_ = new bit_vector_stat();
                break;
            case Config::FAST:
                W_ = new wavelet_tree_fast(bits_per_char_W_);
                last_ = new bit_vector_stat();
                break;
            case Config::SMALL:
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

        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load the BOSS table." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

//
//
// HELPER QUERY FUNCTIONS
//
//

/**
 * Uses the object's array W, a given position i in W and a character c
 * from the alphabet and returns the number of occurences of c in W up to
 * position i.
 */
uint64_t BOSS::rank_W(uint64_t i, TAlphabet c) const {
    assert(i < W_->size());

    return i == 0 ? 0 : W_->rank(c, i) - (c == 0);
}

/**
 * Uses the array W and gets a count i and a character c from
 * the alphabet and returns the position of the i-th occurence of c in W.
 */
uint64_t BOSS::select_W(uint64_t i, TAlphabet c) const {
    assert(i + (c == 0) <= W_->rank(c, W_->size() - 1));

    return i == 0 ? 0 : W_->select(c, i + (c == 0));
}

// get prev character without optimizations (via rank/select calls)
inline uint64_t get_prev(const wavelet_tree &W, uint64_t i, TAlphabet c) {
    assert(i);

    if (W[i] == c)
        return i;

    uint64_t r = W.rank(c, i);
    return r ? W.select(c, r) : 0;
}

/**
 * For characters |first| and |second|, return the last occurrence
 * of them in W[1..i], i.e. max(pred_W(i, first), pred_W(i, second)).
 */
uint64_t BOSS::pred_W(uint64_t i, TAlphabet first, TAlphabet second) const {
    CHECK_INDEX(i);

    if (first == second) {
        uint64_t prev = W_->prev(i, first);
        return prev < W_->size() ? prev : 0;
    }

    // trying to avoid calls of succ_W
    uint64_t max_iter;
    if (dynamic_cast<const wavelet_tree_stat*>(W_)
            || dynamic_cast<const wavelet_tree_fast*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_STAT;
    } else if (dynamic_cast<const wavelet_tree_dyn*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_DYN;
    } else if (dynamic_cast<const wavelet_tree_small*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_SMALL;
    } else {
        assert(false);
        max_iter = 0;
    }

    uint64_t end = i > max_iter
                    ? i - max_iter
                    : 0;

    while (i > end) {
        if (get_W(i) == first)
            return i;
        if (get_W(i) == second)
            return i;
        i--;
    }

    if (!i)
        return 0;

    uint64_t select_first = get_prev(*W_, i, first);
    uint64_t select_second = get_prev(*W_, i, second);

    if (select_second == select_first) {
        return 0;
    } else if (select_first > select_second) {
        return select_first;
    } else {
        return select_second;
    }
}

/**
 * Return the position of the first occurrence of |c| in W[i..N].
 */
uint64_t BOSS::succ_W(uint64_t i, TAlphabet c) const {
    CHECK_INDEX(i);

    return W_->next(i, c);
}

// get next character without optimizations (via rank/select calls)
inline uint64_t get_next(const wavelet_tree &W, uint64_t i, TAlphabet c) {
    assert(i);

    uint64_t r = W.rank(c, i - 1) + 1;
    if (r <= W.rank(c, W.size() - 1)) {
        return W.select(c, r);
    } else {
        return W.size();
    }
}

/**
 * For characters |first| and |second|, return the first occurrence
 * of them in W[i..N], i.e. min(succ_W(i, first), succ_W(i, second)).
 */
std::pair<uint64_t, TAlphabet>
BOSS::succ_W(uint64_t i, TAlphabet first, TAlphabet second) const {
    CHECK_INDEX(i);

    if (first == second) {
        auto next = succ_W(i, first);
        return std::make_pair(next, next < W_->size() ? first : 0);
    }

    // trying to avoid calls of succ_W
    uint64_t max_iter;
    if (dynamic_cast<const wavelet_tree_stat*>(W_)
            || dynamic_cast<const wavelet_tree_fast*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_STAT;
    } else if (dynamic_cast<const wavelet_tree_dyn*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_DYN;
    } else if (dynamic_cast<const wavelet_tree_small*>(W_)) {
        max_iter = MAX_ITER_WAVELET_TREE_SMALL;
    } else {
        assert(false);
        max_iter = 0;
    }

    uint64_t end = std::min(W_->size(), i + max_iter);

    while (i < end) {
        if (get_W(i) == first)
            return std::make_pair(i, first);
        if (get_W(i) == second)
            return std::make_pair(i, second);
        i++;
    }

    if (i == W_->size())
        return std::make_pair(W_->size(), 0);

    uint64_t select_first = get_next(*W_, i, first);
    uint64_t select_second = get_next(*W_, i, second);

    if (select_second == select_first) {
        return std::make_pair(W_->size(), 0);
    } else if (select_first < select_second) {
        return std::make_pair(select_first, first);
    } else {
        return std::make_pair(select_second, second);
    }
}

/**
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
uint64_t BOSS::rank_last(uint64_t i) const {
    assert(i < last_->size());

    return i == 0 ? 0 : last_->rank1(i);
}

/**
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
uint64_t BOSS::select_last(uint64_t i) const {
    assert(i <= last_->num_set_bits());

    return i == 0 ? 0 : last_->select1(i);
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
uint64_t BOSS::pred_last(uint64_t i) const {
    if (!i)
        return 0;

    CHECK_INDEX(i);

    uint64_t prev = last_->prev1(i);

    return prev < last_->size() ? prev : 0;
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
uint64_t BOSS::succ_last(uint64_t i) const {
    CHECK_INDEX(i);

    return last_->next1(i);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character.
 */
uint64_t BOSS::bwd(uint64_t i) const {
    CHECK_INDEX(i);

    if (i == 1)
        return 1;

    uint64_t node_rank = rank_last(i - 1) + 1;

    // get value of last position in node i
    TAlphabet c = get_node_last_value(i);
    // get the offset for the last position in node i
    uint64_t offset = F_[c];
    // compute the offset for this position in W and select it
    return select_W(node_rank - rank_last(offset), c);
}

/**
 * This functions gets a position i reflecting the r-th occurence of the corresponding
 * character c in W and returns the position of the r-th occurence of c in last.
 */
uint64_t BOSS::fwd(uint64_t i) const {
    CHECK_INDEX(i);

    // get value of W at position i
    TAlphabet c = get_W(i) % alph_size;
    assert(i == 1 || c != kSentinelCode);
    // get the offset for position c
    uint64_t o = F_[c];
    // get the rank of c in W at position i
    uint64_t r = rank_W(i, c);
    // select the index of the position in last that is rank many positions after offset
    return select_last(rank_last(o) + r);
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
    return std::make_pair(get_node_last_value(i), bwd(i));
}

/**
 * Given an edge index |i| and a character |c|, get the index of the edge with
 * label c outgoing from the same source node if such exists and npos otherwise.
 */
edge_index BOSS::pick_edge(edge_index edge, TAlphabet c) const {
    CHECK_INDEX(edge);
    assert(get_last(edge) && "must be the last outgoing edge");
    assert(c <= alph_size);

    if (c == alph_size)
        return npos;

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
    uint64_t y = succ_W(x + 1, d);

    // iterate over the rest of the incoming edges
    while (x + 1 < y) {
        x = succ_W(x + 1, d + alph_size);
        if (x < y && get_minus_k_value(x, k_ - 1).first == c) {
            return x;
        }
    }
    return npos;
}

void BOSS::call_incoming_to_target(edge_index edge,
                                   std::function<void(edge_index)> callback) const {
    CHECK_INDEX(edge);
    assert(get_W(edge) < alph_size && "must be the first incoming edge");

    callback(edge);

    const TAlphabet d = get_W(edge);

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
 * Given an edge index i, this function returns true if that is
 * the only edge incoming to its target node.
 */
bool BOSS::is_single_incoming(edge_index i) const {
    CHECK_INDEX(i);

    TAlphabet c = get_W(i);

    assert(c != alph_size);

    if (c > alph_size)
        return false;

    // start from the next edge
    i++;

    return i == W_->size()
            || succ_W(i, c, c + alph_size).second != c + alph_size;
}

/**
 * Given an edge index i (first incoming), this function returns
 * the number of edges incoming to its target node.
 */
size_t BOSS::num_incoming_to_target(edge_index x) const {
    CHECK_INDEX(x);

    assert(get_W(x) < alph_size && "must be the first incoming edge");

    if (x + 1 == W_->size())
        return 1;

    if (dynamic_cast<const wavelet_tree_dyn*>(W_)) {
        TAlphabet d = get_W(x);
        uint64_t y = succ_W(x + 1, d);
        return 1 + rank_W(y - 1, d + alph_size) - rank_W(x - 1, d + alph_size);

    } else {
        size_t indeg = 0;
        call_incoming_to_target(x, [&indeg](auto) { indeg++; });
        assert(indeg && "there is always at least one incoming edge");
        return indeg;
    }
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

    uint64_t last_ = *kmer_it + 1 < alph_size
                     ? F_.at(*kmer_it + 1)
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

        uint64_t last_target = pred_W(last_, s, s + alph_size);
        if (last_target > 0) {
            if (rank_last(last_target - 1) < rank_last(last_ - 1))
                shift = 0;
            last_ = fwd(last_target);
            continue;
        }
        assert(s > 0);

        last_target = succ_W(last_, s, s + alph_size).first;

        if (last_target < W_->size()) {
            last_ = fwd(last_target);
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
 * source node.
 */
std::vector<TAlphabet> BOSS::get_node_seq(edge_index k_node) const {
    CHECK_INDEX(k_node);

    std::vector<TAlphabet> ret(k_, get_node_last_value(k_node));

    for (int curr_k = k_ - 2; curr_k >= 0; --curr_k) {
        CHECK_INDEX(k_node);

        k_node = bwd(k_node);
        ret[curr_k] = get_node_last_value(k_node);
    }

    return ret;
}

/**
 * Given a node index k_node, this function returns the k-mer sequence of the
 * node as a string.
 */
std::string BOSS::get_node_str(edge_index k_node) const {
    CHECK_INDEX(k_node);
    return decode(get_node_seq(k_node));
}

void BOSS::map_to_edges(const std::string &sequence,
                        const std::function<void(edge_index)> &callback,
                        const std::function<bool()> &terminate,
                        const std::function<bool()> &skip) const {
    map_to_edges(encode(sequence), callback, terminate, skip);
}

std::vector<edge_index> BOSS::map_to_edges(const std::string &sequence) const {
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

            edge = fwd(edge);
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

bool BOSS::find(const std::string &sequence,
                double kmer_discovery_fraction) const {
    size_t kmer_size = k_ + 1;

    if (sequence.length() < kmer_size)
        return false;

    const size_t num_kmers = sequence.length() - kmer_size + 1;
    const size_t max_kmers_missing = num_kmers * (1 - kmer_discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    map_to_edges(sequence,
        [&](edge_index edge) {
            if (edge) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; }
    );
    return num_kmers_missing <= max_kmers_missing;
}

bool BOSS::find(const std::string &sequence,
                double kmer_discovery_fraction,
                size_t kmer_mapping_mode) const {
    if (!kmer_mapping_mode)
        return find(sequence, kmer_discovery_fraction);

    size_t kmer_size = k_ + 1;

    if (sequence.length() < kmer_size)
        return false;

    const size_t num_kmers = sequence.length() - kmer_size + 1;
    const size_t max_kmers_missing = num_kmers * (1 - kmer_discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    auto seq_encoded = encode(sequence);

    std::vector<size_t> skipped_kmers;
    skipped_kmers.reserve(seq_encoded.size());

    for (size_t i = 0; i < seq_encoded.size() - kmer_size + 1; ++i) {
        auto edge = map_to_edge(seq_encoded.data() + i,
                                seq_encoded.data() + i + kmer_size);
        if (edge) {
            num_kmers_discovered++;
        } else {
            num_kmers_missing++;
        }

        if (num_kmers_missing > max_kmers_missing
                || num_kmers_discovered >= min_kmers_discovered)
            return num_kmers_discovered >= min_kmers_discovered;

        while (edge && i + kmer_size < seq_encoded.size()) {
            edge = fwd(edge);
            edge = pick_edge(edge, seq_encoded[i + kmer_size]);
            if (edge) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }

            if (num_kmers_missing > max_kmers_missing
                    || num_kmers_discovered >= min_kmers_discovered)
                return num_kmers_discovered >= min_kmers_discovered;

            i++;
        }

        if (kmer_discovery_fraction < 1
            && (kmer_mapping_mode == 1
                || kmer_size * (max_kmers_missing - num_kmers_missing)
                    > min_kmers_discovered - num_kmers_discovered)) {
            size_t i_old = i;

            // jump over the missing k-mer
            i += kmer_size - 1;

            // Save skipped kmers for the end
            while (++i_old <= i && i_old < seq_encoded.size() - kmer_size + 1) {
                skipped_kmers.push_back(i_old);
            }
        }
    }

    for (size_t j = 0; j < skipped_kmers.size(); ++j) {
        size_t i = skipped_kmers[j];

        auto edge = map_to_edge(seq_encoded.data() + i,
                                seq_encoded.data() + i + kmer_size);
        if (edge) {
            num_kmers_discovered++;
        } else {
            num_kmers_missing++;
        }

        if (num_kmers_missing > max_kmers_missing
                || num_kmers_discovered >= min_kmers_discovered)
            return num_kmers_discovered >= min_kmers_discovered;

        while (edge && j + 1 < skipped_kmers.size()
                    && i + 1 == skipped_kmers[j + 1]) {
            edge = fwd(edge);
            edge = pick_edge(edge, seq_encoded[i + kmer_size]);
            if (edge) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }

            if (num_kmers_missing > max_kmers_missing
                    || num_kmers_discovered >= min_kmers_discovered)
                return num_kmers_discovered >= min_kmers_discovered;

            i++;
            j++;
        }
    }

    return num_kmers_missing <= max_kmers_missing;
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
 * This function gets a value of the alphabet c and updates the offset of
 * all following values by +1 is positive is true and by -1 otherwise.
 */
void BOSS::update_F(TAlphabet c, int value) {
    assert(c < alph_size);
    assert(std::abs(value) == 1);

    for (TAlphabet i = c + 1; i < alph_size; i++) {
        F_[i] += value;
    }
}

TAlphabet BOSS::encode(char s) const {
    assert(kmer_extractor_.encode(kSentinel) != kSentinelCode);
    assert(kmer_extractor_.encode(s) <= alph_size);
    return kmer_extractor_.encode(s);
}

std::vector<TAlphabet> BOSS::encode(const std::string &sequence) const {
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

void BOSS::switch_state(Config::StateType new_state) {

    //std::cerr << "switching state from " << this->state << " to " << state << std::endl;
    if (state == new_state)
        return;

    switch (new_state) {
        case Config::STAT: {
            convert<wavelet_tree_stat, bit_vector_stat>(&W_, &last_);
            break;
        }
        case Config::SMALL: {
            convert<wavelet_tree_small, bit_vector_small>(&W_, &last_);
            break;
        }
        case Config::FAST: {
            convert<wavelet_tree_fast, bit_vector_stat>(&W_, &last_);
            break;
        }
        case Config::DYN: {
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
    for (uint64_t i = 0; i < W_->size(); ++i) {
        os << " " << static_cast<uint64_t>(get_W(i));
    }
    os << std::endl;
}

void BOSS::print(std::ostream &os) const {
    assert(is_valid());
    auto vertex_header = std::string("Vertex");
    vertex_header.resize(k_, ' ');

    os << "Index" << "\t" << "L"
                  << "\t" << vertex_header
                  << "\t" << "W" << std::endl;

    for (uint64_t i = 1; i < W_->size(); i++) {
        assert(get_W(i) != alph_size);
        os << i << "\t" << get_last(i)
                << "\t" << get_node_str(i)
                << "\t" << decode(get_W(i) % alph_size)
                        << (get_W(i) > alph_size
                                ? "-"
                                : "")
                        << std::endl;
    }
}

void BOSS::print_adj_list(std::ostream &os) const {
    for (uint64_t edge = 1; edge < W_->size(); ++edge) {
        os << 1 + rank_last(fwd(edge) - 1)
           << " ";
        if (get_last(edge))
            os << "\n";
    }
}

///////////////
// Construct //
///////////////

// add a full sequence to the graph
void BOSS::add_sequence(const std::string &seq,
                        bool try_extend,
                        std::vector<uint64_t> *edges_inserted) {
    if (seq.size() < k_ + 1)
        return;

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

            uint64_t source;

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
                            std::vector<uint64_t> *edges_inserted) {
    CHECK_INDEX(source_node);
    assert(source_node_kmer);
    assert(std::vector<TAlphabet>(source_node_kmer, source_node_kmer + k_)
                                                == get_node_seq(source_node));
    assert(c < alph_size);

    // get range of identical nodes (without W) pos current end position
    uint64_t begin = pred_last(source_node - 1) + 1;
    uint64_t end = succ_last(source_node) + 1;

    // get position of the first occurence of c or c- in W after p
    uint64_t prev_c_pos = pred_W(end - 1, c, c + alph_size);
    // if the edge already exists, traverse it and return the index
    if (prev_c_pos >= begin)
        return fwd(prev_c_pos);

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
        uint64_t inserted = insert_edge(c + alph_size, begin, end);
        if (edges_inserted && inserted)
            edges_inserted->push_back(inserted);

        return fwd(prev_c_pos);
    }

    // The new edge will be the first incoming for its target node,
    // and therefore the new edge will be marked by c (not c-)

    // adding a new node can influence one of the following nodes sharing the k-1 suffix
    // get position of the first occurence of c after p (including p + 1)
    uint64_t first_c = end < W_->size()
                       ? succ_W(end, c)
                       : W_->size();

    bool the_only_incoming = true;
    if (first_c < W_->size()
            && !(the_only_incoming = !compare_node_suffix(first_c, source_node_kmer))) {
        // The inserted edge will not be the only incoming for its target node.
        // Relabel the next incoming edge from c to c- since
        // the new edge is the first incoming for that target node.
        W_->set(first_c, c + alph_size);
    }

    // insert the edge
    uint64_t inserted = insert_edge(c, begin, end);
    if (edges_inserted && inserted)
        edges_inserted->push_back(inserted);

    // Add sentinel if the target node is the new dead-end
    if (!the_only_incoming)
        return fwd(first_c + (inserted > 0));

    uint64_t sentinel_pos = select_last(rank_last(F_[c]) + rank_W(begin - 1, c)) + 1;

    update_F(c, +1);
    W_->insert(sentinel_pos, kSentinelCode);
    last_->insert_bit(sentinel_pos, true);

    if (edges_inserted)
        edges_inserted->push_back(sentinel_pos);

    assert((*W_)[0] == 0);

    return sentinel_pos;
}


uint64_t BOSS::insert_edge(TAlphabet c, uint64_t begin, uint64_t end) {
    assert(c != alph_size);
    assert(c < 2 * alph_size);

    if (begin > 1 && get_W(begin) == kSentinelCode) {
        // the source node is the dead-end with outgoing sentinel
        // replace this sentinel with the proper label
        W_->set(begin, c);
        return 0;
    } else {
        // the source node already has some outgoing edges

        // find the exact position of the new edge
        uint64_t pos = begin;
        while (pos < end && get_W(pos) % alph_size < c % alph_size) {
            pos++;
        }

        // insert the new edge
        update_F(get_node_last_value(begin), +1);
        last_->insert_bit(begin, false);
        W_->insert(pos, c);

        assert(pos);
        return pos;
    }
}


// Given an edge list, remove them from the BOSS graph.
// TODO: fix the implementation (anchoring the isolated nodes)
void BOSS::erase_edges_dyn(const std::set<edge_index> &edges) {
    uint64_t shift = 0;

    for (edge_index edge : edges) {
        assert(edge >= shift);
        uint64_t edge_id = edge - shift;

        uint64_t d = get_W(edge_id);
        if (d < alph_size && edge_id + 1 < W_->size()) {
            //fix W array
            auto [next, d_next] = succ_W(edge_id + 1, d, d + alph_size);
            if (d_next == d + alph_size)
                W_->set(next, d);
        }
        W_->remove(edge_id);
        update_F(get_node_last_value(edge_id), -1);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        if (get_last(edge_id) && (edge >= shift + 1)
                              && !get_last(edge_id - 1)) {
            last_->delete_bit(edge_id - 1);
        } else {
            last_->delete_bit(edge_id);
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
    sdsl::bit_vector new_last(last_->size() - num_edges_to_remove, false);

    for (uint64_t i = 0, new_i = 0; i < edges_to_remove_mask.size(); ++i) {
        if (!edges_to_remove_mask[i]) {
            new_last[new_i++] = get_last(i);
        } else {
            if (get_last(i) && new_i > 1 && !new_last[new_i - 1])
                new_last[new_i - 1] = 1;
        }
    }
    delete last_;
    last_ = new bit_vector_stat(std::move(new_last));

    // update W
    sdsl::int_vector<> new_W(W_->size() - num_edges_to_remove, 0, bits_per_char_W_);
    sdsl::bit_vector first_removed(alph_size, false);

    for (uint64_t i = 0, new_i = 0; i < edges_to_remove_mask.size(); ++i) {
        TAlphabet c = get_W(i);
        if (edges_to_remove_mask[i]) {
            if (c < alph_size)
                first_removed[c] = true;
        } else {
            assert(c != alph_size);
            if (c > alph_size && first_removed[c % alph_size]) {
                new_W[new_i++] = c % alph_size;
            } else {
                new_W[new_i++] = c;
            }
            first_removed[c % alph_size] = false;
        }
    }
    delete W_;
    W_ = new wavelet_tree_stat(bits_per_char_W_, std::move(new_W));

    // update F
    TAlphabet c = 0;
    uint64_t count = 0;
    for (uint64_t i = 1; i <= F_.back(); ++i) {
        while (i > F_[c] && c < alph_size) {
            F_[c++] = count;
        }
        if (!edges_to_remove_mask[i])
            count++;
    }
    while (c < alph_size) {
        F_[c++] = count;
    }

    state = Config::STAT;

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
            path.push_back(pred_last(fwd(path.back()) - 1) + 1);
            pre_visit(path.back());
        }

        // traverse the path backwards to the next branching node
        while (path.size() > 1 && get_last(path.back())) {
            post_visit(path.back());
            path.pop_back();
        }

        // explore the next edge outgoing from the current branching node
        if (path.size() > 1) {
            post_visit(path.back());
            path.back()++;
            pre_visit(path.back());
        }
    } while (path.size() > 1);

    post_visit(path.back());
}

/**
 * Traverse the entire dummy subtree
 * and find all redundant dummy edges.
 */
template <typename Array, typename U>
void traverse_dummy_edges(const BOSS &graph,
                          edge_index subtree_root,
                          size_t check_depth,
                          Array *redundant_mask,
                          Array *traversed_mask,
                          U *num_dummy_traversed,
                          bool verbose) {
    assert(!redundant_mask || redundant_mask->size() == graph.get_W().size());
    assert(!traversed_mask || traversed_mask->size() == graph.get_W().size());

    if (!check_depth)
        return;

    std::vector<bool> redundant_path;
    uint64_t num_edges_traversed = 0;

    // start traversal in the given node
    graph.edge_DFT(subtree_root,
        [&](edge_index) {
            redundant_path.push_back(true);
            num_edges_traversed++;
        },
        [&](edge_index edge) {
            assert(graph.get_W(edge) < graph.alph_size);

            if (traversed_mask)
                (*traversed_mask)[edge] = 1;

            if (redundant_path.back())
                (*redundant_mask)[edge] = 2;

            redundant_path.pop_back();
        },
        [&](edge_index edge) {
            if (redundant_path.size() == check_depth) {
                if (!redundant_mask
                        || (!(*redundant_mask)[edge]
                                && graph.is_single_incoming(edge))) {
                    // the last dummy edge is not redundant and hence the
                    // entire path has to remain in the graph
                    redundant_path.assign(redundant_path.size(), false);
                }
                return true;
            } else {
                assert(graph.is_single_incoming(edge));
                return false;
            }
        }
    );
    *num_dummy_traversed += num_edges_traversed;
    if (verbose) {
        std::cout << std::string("Source dummy edges traversed: ")
                        + std::to_string(*num_dummy_traversed) + "\n" << std::flush;
    }
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
                              bool verbose) {
    assert(!redundant_mask || redundant_mask->size() == graph.get_W().size());
    assert(!traversed_mask || traversed_mask->size() == graph.get_W().size());

    if (traversed_mask)
        (*traversed_mask)[1] = true;

    if (graph.get_last(1))
        return 1;

    if (num_threads <= 1 || graph.get_k() <= 1) {
        // assume the main dummy source node already traversed
        uint64_t num_dummy_traversed = 1;
        // start traversal in the main dummy source node
        uint64_t root = 2;
        do {
            traverse_dummy_edges(graph, root, graph.get_k(),
                                 redundant_mask,
                                 traversed_mask,
                                 &num_dummy_traversed,
                                 verbose);
        } while (!graph.get_last(root++));

        return num_dummy_traversed;
    }

    ThreadPool pool(num_threads);

    std::unique_ptr<std::vector<char>> edges_threadsafe;
    if (redundant_mask || traversed_mask)
        edges_threadsafe.reset(new std::vector<char>(graph.get_W().size(), 0));

    // assume the main dummy source node already traversed
    std::atomic<uint64_t> num_dummy_traversed { 1 };

    const size_t tree_split_depth = std::min(size_t(6), graph.get_k() / 2);

    // run traversal for subtree at depth |tree_split_depth| in parallel
    uint64_t root = 2;
    size_t depth = 0;
    do {
        graph.edge_DFT(root,
            [&](edge_index) { depth++; },
            [&](edge_index) { depth--; },
            [&](edge_index edge) {
                if (depth == tree_split_depth) {
                    pool.enqueue([](const auto& ...args) {
                                    traverse_dummy_edges(args...);
                                 },
                                 std::ref(graph),
                                 edge, 1 + graph.get_k() - tree_split_depth,
                                 edges_threadsafe.get(),
                                 edges_threadsafe.get(),
                                 &num_dummy_traversed, verbose);
                    num_dummy_traversed--;
                    return true;
                } else {
                    return false;
                }
            }
        );
        assert(!depth);
    } while (!graph.get_last(root++));

    pool.join();
    if (verbose) {
        std::cout << "All subtrees are done" << std::endl;
    }

    if (edges_threadsafe.get()) {
        for (uint64_t i = 0; i < edges_threadsafe->size(); ++i) {
            if ((*edges_threadsafe)[i] && traversed_mask)
                (*traversed_mask)[i] = true;
            if ((*edges_threadsafe)[i] == 2 && redundant_mask)
                (*redundant_mask)[i] = true;
        }
        edges_threadsafe.reset();
    }

    // propogate the results from the subtrees
    root = 2;
    do {
        traverse_dummy_edges(graph, root, tree_split_depth,
                             redundant_mask,
                             traversed_mask,
                             &num_dummy_traversed,
                             verbose);
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
    sdsl::bit_vector redundant_dummy_edges_mask(W_->size(), false);

    if (source_dummy_edges) {
        (*source_dummy_edges) = sdsl::bit_vector();
        (*source_dummy_edges) = sdsl::bit_vector(W_->size(), false);
        (*source_dummy_edges)[1] = true;
    }

    if (get_last(1))
        return redundant_dummy_edges_mask;

    switch_state(Config::STAT);

    auto num_dummy_traversed = traverse_dummy_edges(
        *this, &redundant_dummy_edges_mask, source_dummy_edges, num_threads, verbose
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

    return redundant_dummy_edges_mask;
}

uint64_t BOSS::mark_source_dummy_edges(sdsl::bit_vector *mask,
                                       size_t num_threads,
                                       bool verbose) const {
    assert(!mask || mask->size() == W_->size());

    return traverse_dummy_edges(*this, NULL, mask, num_threads, verbose);
}

uint64_t BOSS::mark_sink_dummy_edges(sdsl::bit_vector *mask) const {
    if (!mask)
        return rank_W(num_edges(), 0) - 1;

    assert(mask->size() == W_->size());

    uint64_t num_dummy_sink_edges = 0;

    // skip the main dummy source
    for (uint64_t i = 2; i < W_->size(); ++i) {
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
void BOSS::merge(const BOSS &other) {
    other.call_sequences([&](const std::string &sequence, auto&&) {
        add_sequence(sequence, true);
    });
}

void BOSS::call_start_edges(Call<edge_index> callback) const {
    // start traversal in the main dummy source node and traverse the tree
    // of source dummy k-mers to depth k + 1

    // check if the dummy tree is not empty
    if (is_single_outgoing(1))
        return;

    // run traversal for subtree
    uint64_t root = 2;
    size_t depth = 0;
    do {
        edge_DFT(root,
            [&](edge_index) { depth++; },
            [&](edge_index) { depth--; },
            [&](edge_index edge) {
                if (depth < k_)
                    return false;

                if (depth == k_)
                    return !is_single_incoming(edge);

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
// If multiple outgoing edges are found, set |*i| to the first and return false.
bool masked_pick_single_outgoing(const BOSS &boss,
                                 uint64_t *i,
                                 const bitmap *subgraph_mask) {
    assert(i && *i);
    assert(boss.get_last(*i));
    assert(!subgraph_mask || subgraph_mask->size() == boss.num_edges() + 1);

    // in boss, at least one outgoing edge always exists
    if (!subgraph_mask)
        return boss.is_single_outgoing(*i);

    bool edge_detected = false;
    uint64_t j = *i;
    do {
        if ((*subgraph_mask)[j]) {
            // there are multiple outgoing edges
            if (edge_detected)
                return false;

            // this is the first outgoing edge
            edge_detected = true;
            *i = j;
        }
    } while (--j > 0 && !boss.get_last(j));

    // return true of there is exactly one outgoing edge
    if (edge_detected)
        return true;

    // no outgoing
    *i = 0;
    return false;
}

// If a single incoming edge is found, write it to |*i| and return true.
// If no incoming edges are found, set |*i| to 0 and return false.
// If multiple incoming edges are found, set |*i| to the first and return false.
bool masked_pick_single_incoming(const BOSS &boss,
                                 uint64_t *i,
                                 const bitmap *subgraph_mask) {
    assert(i && *i);
    assert(boss.get_W(*i) < boss.alph_size);
    assert(!subgraph_mask || subgraph_mask->size() == boss.num_edges() + 1);

    // in boss, at least one incoming edge always exists
    if (!subgraph_mask)
        return boss.is_single_incoming(*i);

    auto d = boss.get_W(*i);
    auto j = *i;

    TAlphabet d_next;
    bool edge_detected = false;
    do {
        if ((*subgraph_mask)[j]) {
            // there are multiple incoming edges
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

// traverse graph from the specified (k+1)-mer/edge and call
// all paths reachable from it
void call_paths(const BOSS &boss,
                edge_index starting_kmer,
                BOSS::Call<std::vector<edge_index>&&,
                           std::vector<TAlphabet>&&> callback,
                bool split_to_unitigs,
                sdsl::bit_vector *discovered_ptr,
                sdsl::bit_vector *visited_ptr,
                ProgressBar &progress_bar,
                const bitmap *subgraph_mask);

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void BOSS::call_paths(Call<std::vector<edge_index>&&,
                           std::vector<TAlphabet>&&> callback,
                      bool split_to_unitigs,
                      const bitmap *subgraph_mask) const {
    assert(!subgraph_mask || subgraph_mask->size() == W_->size());

    // keep track of reached edges
    sdsl::bit_vector discovered(W_->size(), false);
    if (subgraph_mask) {
        subgraph_mask->add_to(&discovered);
        discovered.flip();
    }
    discovered[0] = true;
    // keep track of edges that are already included in covering paths
    sdsl::bit_vector visited = discovered;

    ProgressBar progress_bar(visited.size() - sdsl::util::cnt_one_bits(visited),
                             "Traverse BOSS",
                             std::cerr, !utils::get_verbose());

    // process source dummy edges first
    //
    //  .____
    //
    if (!subgraph_mask) {
        for (uint64_t i = succ_last(1); i >= 1; --i) {
            if (!visited[i])
                ::call_paths(*this, i, callback, split_to_unitigs,
                             &discovered, &visited, progress_bar, nullptr);
        }

    } else {
        call_zeros(visited, [&](uint64_t i) {
            if (!get_last(i))
                return;

            // check if |i| has incoming edges
            auto j = bwd(i);
            // check if =1 or >1
            if (masked_pick_single_incoming(*this, &j, subgraph_mask) || j)
                return;

            do {
                if (!visited[i])
                    ::call_paths(*this, i, callback, split_to_unitigs,
                                 &discovered, &visited, progress_bar, subgraph_mask);
            } while (--i > 0 && !get_last(i));
        });
    }

    // then all forks
    //  ____.____
    //       \___
    //
    call_zeros(visited, [&](uint64_t i) {
        if (!get_last(i))
            return;

        if (masked_pick_single_outgoing(*this, &i, subgraph_mask) || !i)
            return;

        do {
            if (!visited[i])
                ::call_paths(*this, i, callback, split_to_unitigs,
                             &discovered, &visited, progress_bar, subgraph_mask);
        } while (--i > 0 && !get_last(i));
    });

    // process all the cycles left that have not been traversed
    call_zeros(visited, [&](uint64_t i) {
        ::call_paths(*this, i, callback, split_to_unitigs,
                     &discovered, &visited, progress_bar, subgraph_mask);
    });
}

struct Edge {
    BOSS::edge_index id;
    std::vector<TAlphabet> source_kmer;
};

void call_paths(const BOSS &boss,
                edge_index starting_kmer,
                BOSS::Call<std::vector<edge_index>&&,
                           std::vector<TAlphabet>&&> callback,
                bool split_to_unitigs,
                sdsl::bit_vector *discovered_ptr,
                sdsl::bit_vector *visited_ptr,
                ProgressBar &progress_bar,
                const bitmap *subgraph_mask) {
    assert(discovered_ptr && visited_ptr);

    auto &discovered = *discovered_ptr;
    auto &visited = *visited_ptr;
    // store all branch nodes on the way
    std::vector<TAlphabet> kmer;
    discovered[starting_kmer] = true;
    std::deque<Edge> edges { { starting_kmer, boss.get_node_seq(starting_kmer) } };

    // keep traversing until we have worked off all branches from the queue
    while (!edges.empty()) {
        std::vector<uint64_t> path;
        uint64_t edge = edges.back().id;
        auto sequence = std::move(edges.back().source_kmer);
        edges.pop_back();

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (!visited[edge]) {
            assert(edge > 0 && discovered[edge]);
            assert(!subgraph_mask || (*subgraph_mask)[edge]);

            // visit the edge
            sequence.push_back(boss.get_W(edge) % boss.alph_size);
            path.push_back(edge);
            visited[edge] = true;
            ++progress_bar;

            // stop traversing if the next node is a dummy sink
            if (!sequence.back())
                break;

            // stop traversing if we call unitigs and this
            // is not the only incoming edge
            auto d = boss.get_W(edge) % boss.alph_size;
            auto j = boss.pred_W(edge, d, d);
            bool continue_traversal = !split_to_unitigs
                || masked_pick_single_incoming(boss, &j, subgraph_mask);
            assert(j);

            // make one traversal step
            assert(boss.get_W(edge));
            edge = boss.fwd(edge);

            // traverse if there is only one outgoing edge
            auto is_single = masked_pick_single_outgoing(boss, &edge, subgraph_mask);

            if (continue_traversal && is_single) {
                discovered[edge] = true;
                continue;
            } else if (!edge) {
                break;
            }

            kmer.assign(sequence.end() - boss.get_k(), sequence.end());
            edge_index next_edge = 0;

            // loop over outgoing edges
            do {
                if (!next_edge && !split_to_unitigs && !visited[edge]) {
                    // save the edge for visiting if we extract arbitrary paths
                    discovered[edge] = true;
                    next_edge = edge;
                } else if (!discovered[edge]) {
                    // discover other edges
                    discovered[edge] = true;
                    edges.push_back({ edge, kmer });
                }
            } while (--edge > 0 && !boss.get_last(edge));

            // stop traversing this sequence if the next edge was not selected
            if (!next_edge)
                break;

            // pick the last outgoing but not yet visited
            // edge and continue traversing the graph
            edge = next_edge;
        }

        if (path.size())
            callback(std::move(path), std::move(sequence));
    }
}

void BOSS::call_sequences(Call<std::string&&, std::vector<uint64_t>&&> callback,
                          const bitmap *subgraph_mask) const {

    call_paths([&](auto&& edges, auto&& path) {
        assert(path.size());

        auto begin = std::find_if(path.begin(), path.end(),
                                  [&](auto c) { return c != kSentinelCode; });
        auto end = path.end() - (path.back() == kSentinelCode);

        if (begin + k_ + 1 > end)
            return;

        assert(std::all_of(begin, end, [](TAlphabet c) { return c != kSentinelCode; }));

        std::string sequence(end - begin, '\0');
        std::transform(begin, end, sequence.begin(),
                       [&](auto c) { return BOSS::decode(c % alph_size); });

        std::copy(edges.begin() + (begin - path.begin()),
                  edges.end() - (path.end() - end),
                  edges.begin());
        edges.resize(edges.size()
                        - (begin - path.begin())
                        - (path.end() - end));

        assert(end - begin - k_ == edges.size());

        callback(std::move(sequence), std::move(edges));

    }, false, subgraph_mask);
}

void BOSS::call_unitigs(Call<std::string&&, std::vector<edge_index>&&> callback,
                        size_t min_tip_size,
                        const bitmap *subgraph_mask) const {
    call_paths([&](auto&& edges, auto&& path) {
        assert(path.size());

        auto begin = std::find_if(path.begin(), path.end(),
                                  [&](auto c) { return c != kSentinelCode; });
        auto end = path.end() - (path.back() == kSentinelCode);

        if (begin + k_ + 1 > end)
            return;

        assert(std::all_of(begin, end, [](TAlphabet c) { return c != kSentinelCode; }));

        std::string sequence(end - begin, '\0');
        std::transform(begin, end, sequence.begin(),
                       [&](auto c) { return BOSS::decode(c % alph_size); });

        auto first_edge = edges.front();
        auto last_edge = edges.back();

        std::copy(edges.begin() + (begin - path.begin()),
                  edges.end() - (path.end() - end),
                  edges.begin());
        edges.resize(edges.size()
                        - (begin - path.begin())
                        - (path.end() - end));

        assert(end - begin - k_ == edges.size());

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

        uint64_t last_fwd = 0;

        // if the last node has multiple outgoing edges,
        // it is clearly neither a sink tip nor a source tip.
        if (path.back() != kSentinelCode
                && !masked_pick_single_outgoing(*this,
                                                &(last_fwd = fwd(last_edge)),
                                                subgraph_mask)
                && last_fwd) {
            callback(std::move(sequence), std::move(edges));
            return;
        }

        uint64_t first_bwd = 0;

        // if the first node has multiple incoming edges,
        // it is clearly neither a source tip nor a sink tip.
        // TODO: This doesn't work properly if graph has redundant dummy edges.
        //       Make sure there are no redundant dummy edges when this
        //       function is called.
        if (path.front() != kSentinelCode
                && !masked_pick_single_incoming(*this,
                                                &(first_bwd = bwd(first_edge)),
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

    }, true, subgraph_mask);
}

/**
 * Traverse boss graph and call all its edges
 * except for the dummy source of sink ones
 */
void BOSS::call_kmers(Call<edge_index, const std::string&> callback) const {
    sdsl::bit_vector visited(W_->size(), false);

    // store all branch nodes on the way
    std::queue<std::pair<edge_index, std::string>> branchnodes;

    // start from the second edge (skip dummy main source)
    for (uint64_t i = 2; i < W_->size(); ++i) {
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

                // stop traversing if it's a sink
                if (!get_W(edge))
                    break;

                // traverse if there is only one outgoing edge
                if (is_single_outgoing(edge)) {
                    auto next_edge = fwd(edge);

                    kmer.back() = decode(get_node_last_value(next_edge));
                    if (kmer.front() != BOSS::kSentinel)
                        callback(edge, kmer);

                    edge = next_edge;
                    std::copy(kmer.begin() + 1, kmer.end(), kmer.begin());

                } else {
                    // loop over outgoing edges
                    do {
                        assert(get_W(edge));

                        auto next_edge = fwd(edge);

                        kmer.back() = decode(get_node_last_value(next_edge));
                        if (kmer.front() != BOSS::kSentinel)
                            callback(edge, kmer);

                        if (!visited[next_edge])
                            branchnodes.push({ next_edge, kmer.substr(1) + '\0' });

                    } while (--edge > 1 && !get_last(edge));

                    break;
                }
            }
        }
    }
}

bool BOSS::is_valid() const {
    assert((*W_)[0] == 0);
    assert(W_->size() >= 2);
    assert(get_node_str(1) == std::string(k_, kSentinel) && "First kmer must be dummy");
    assert(get_W(1) == kSentinelCode && "First kmer must be dummy");

    for (uint64_t i = 1; i < W_->size(); i++) {
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
