#include "dbg_succinct.hpp"

/**
 * This class contains a succinct representation of the de bruijn graph
 * following ideas and suggestions presented here:
 * http://link.springer.com/chapter/10.1007/978-3-642-33122-0_18
 *
 * There is also conceptual code available at
 * https://code.google.com/p/csalib/downloads/list
 * that has been used as a reference for this implementation.
 */

#include <cassert>
#include <vector>
#include <stack>
#include <algorithm>
#include <string>
#include <cstdio>

#include <boost/multiprecision/integer.hpp>

#include "dbg_succinct_construct.hpp"
#include "serialization.hpp"
#include "helpers.hpp"

using utils::remove_suffix;
using TAlphabet = DBG_succ::TAlphabet;


#define CHECK_INDEX(idx) \
    assert(idx < W_->size()); \
    assert(idx > 0)
#define CHECK_NODE(idx) \
    assert(idx <= num_nodes()); \
    assert(idx > 0)

const DBG_succ::node_index DBG_succ::npos = 0;

typedef DBG_succ::node_index node_index;
typedef DBG_succ::edge_index edge_index;


DBG_succ::DBG_succ(size_t k)
      : alph_size(kmer_extractor_.alphabet.size()),
        alphabet(kmer_extractor_.alphabet),
        bits_per_char_W_(boost::multiprecision::msb(alph_size - 1) + 2),
        k_(k),
        last_(new bit_vector_dyn()),
        F_(alph_size, 0),
        W_(new wavelet_tree_dyn(bits_per_char_W_)) {

    assert(bits_per_char_W_ <= sizeof(TAlphabet) * 8
            && "Choose type for TAlphabet properly");

    last_->insertBit(0, false);
    W_->insert(0, 0);

    // add the dummy source node
    last_->insertBit(1, true);
    W_->insert(0, 0);
    for (size_t j = 1; j < alph_size; j++) {
        F_[j] = 1;
    }
    assert(is_valid());
}

DBG_succ::DBG_succ(DBGSuccConstructor *builder) : DBG_succ::DBG_succ() {
    assert(builder);

    builder->build_graph(this);
    assert(is_valid());
}

DBG_succ::~DBG_succ() {
    delete W_;
    delete last_;
}

/**
 * Given a pointer to a graph structures G1 and G2, the function compares their elements to the
 * each other. It will perform an element wise comparison of the arrays W, last and
 * F and will only check for identity. If any element differs, the function will return
 * false and true otherwise.
 */
bool DBG_succ::equals_internally(const DBG_succ &other, bool verbose) const {
    // compare size
    if (W_->size() != other.W_->size()) {
        if (verbose)
            std::cout << "sizes of graphs differ"
                      << "\n1: " << W_->size()
                      << "\n2: " << other.W_->size()
                      << std::endl;
        return false;
    }

    assert(F_.size() == other.F_.size());

    // compare last
    for (size_t i = 0; i < W_->size(); ++i) {
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
    for (size_t i = 0; i < W_->size(); ++i) {
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
    for (size_t i = 0; i < F_.size(); ++i) {
        if (get_F(i) != other.get_F(i)) {
            if (verbose)
                std::cout << "F differs at position " << i
                          << "\n1: F[" << i << "] = " << get_F(i)
                          << "\n2: F[" << i << "] = " << other.get_F(i)
                          << std::endl;
            return false;
        }
    }

    return true;
}

/**
 * Check whether graphs store the same data.
 * FYI: this function reconstructs all the kmers, so
 * the complexity is at least O(k x n).
 */
bool DBG_succ::operator==(const DBG_succ &other) const {
    uint64_t i = 1;
    uint64_t j = 1;

    while (i < W_->size() && j < other.W_->size()) {

        std::string first_node, second_node;
        bool first_last, second_last;
        char first_label, second_label;

        do {
            first_node = get_node_str(i);
            first_last = get_last(i);
            first_label = decode(get_W(i));
            i++;
        } while (first_node.find(kSentinel) != std::string::npos && i < W_->size());

        do {
            second_node = other.get_node_str(j);
            second_last = other.get_last(j);
            second_label = decode(other.get_W(j));
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

void DBG_succ::serialize(const std::string &filename) const {
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

void DBG_succ::serialize(std::ofstream &outstream) const {
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

bool DBG_succ::load(const std::string &filename) {
    auto file = remove_suffix(filename, kExtension) + kExtension;

    std::ifstream instream(file, std::ios::binary);

    return load(instream);
}

bool DBG_succ::load(std::ifstream &instream) {
    // if not specified in the file, the default for loading is dynamic
    state = Config::DYN;

    try {
        // load F, k, and state
        F_ = libmaus2::util::NumberSerialisation::deserialiseNumberVector<uint64_t>(instream);
        k_ = load_number(instream);
        state = static_cast<Config::StateType>(load_number(instream));

        if (F_.size() != alph_size)
            return false;

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
            case Config::SMALL:
                W_ = new wavelet_tree_small(bits_per_char_W_);
                last_ = new bit_vector_small();
                break;
        }
        return W_->load(instream) && last_->load(instream);
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load DBG_succ." << std::endl;
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
uint64_t DBG_succ::rank_W(uint64_t i, TAlphabet c) const {
    assert(i < W_->size());

    return i == 0 ? 0 : W_->rank(c, i) - (c == 0);
}

/**
 * Uses the array W and gets a count i and a character c from
 * the alphabet and returns the position of the i-th occurence of c in W.
 */
uint64_t DBG_succ::select_W(uint64_t i, TAlphabet c) const {
    assert(i + (c == 0) <= W_->rank(c, W_->size() - 1));

    return i == 0 ? 0 : W_->select(c, i + (c == 0));
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the last index of a character c preceding in W[1..i].
 */
uint64_t DBG_succ::pred_W(uint64_t i, TAlphabet c) const {
    assert(i < W_->size());

    size_t max_iter = 10;
    if (state == Config::STAT) {
        max_iter = 1000;
    }
    for (size_t t = 0; t < max_iter; ++t, --i) {
        if (i == 0 || get_W(i) == c)
            return i;
    }

    return select_W(rank_W(i, c), c);
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the first index of a character c in W[i..N].
 */
uint64_t DBG_succ::succ_W(uint64_t i, TAlphabet c) const {
    assert(i < W_->size());

    size_t max_iter = 10;
    if (state == Config::STAT) {
        max_iter = 1000;
    }
    for (size_t t = 0; t < max_iter; ++t) {
        if (i + t == W_->size() || get_W(i + t) == c)
            return i + t;
    }

    uint64_t rk = rank_W(i, c);
    if (rk == rank_W(W_->size() - 1, c))
        return W_->size();

    return select_W(rk + 1, c);
}

/**
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
uint64_t DBG_succ::rank_last(uint64_t i) const {
    assert(i < last_->size());

    return i == 0 ? 0 : last_->rank1(i);
}

/**
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
uint64_t DBG_succ::select_last(uint64_t i) const {
    assert(i <= last_->num_set_bits());

    return i == 0 ? 0 : last_->select1(i);
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
uint64_t DBG_succ::pred_last(uint64_t i) const {
    assert(i < last_->size());

    size_t max_iter = 10;
    if (state == Config::STAT) {
        max_iter = 1000;
    }
    for (size_t t = 0; t < max_iter; ++t, --i) {
        if (i == 0 || get_last(i))
            return i;
    }

    return select_last(rank_last(i));
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
uint64_t DBG_succ::succ_last(uint64_t i) const {
    CHECK_INDEX(i);

    size_t max_iter = 10;
    if (state == Config::STAT) {
        max_iter = 1000;
    }
    for (size_t t = 0; t < max_iter; ++t) {
        if (i + t == W_->size() || get_last(i + t))
            return i + t;
    }

    uint64_t next_rank = get_source_node(i);

    assert(next_rank <= last_->num_set_bits());

    return select_last(next_rank);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character.
 */
uint64_t DBG_succ::bwd(uint64_t i) const {
    CHECK_INDEX(i);

    uint64_t node_rank = get_source_node(i);

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
uint64_t DBG_succ::fwd(uint64_t i) const {
    CHECK_INDEX(i);

    // get value of W at position i
    TAlphabet c = get_W(i) % alph_size;
    // get the offset for position c
    uint64_t o = F_[c];
    // get the rank of c in W at position i
    uint64_t r = rank_W(i, c);
    // select the index of the position in last that is rank many positions after offset
    return select_last(rank_last(o) + r);
}

node_index DBG_succ::get_source_node(edge_index i) const {
    CHECK_INDEX(i);
    return rank_last(i - 1) + 1;
}

/**
 * Using the offset structure F this function returns the value of the last
 * position of the source node for edge i.
 */
TAlphabet DBG_succ::get_node_last_value(edge_index i) const {
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
DBG_succ::get_minus_k_value(edge_index i, size_t k) const {
    CHECK_INDEX(i);

    for (; k > 0; --k) {
        i = bwd(i);
    }
    return std::make_pair(get_node_last_value(i), bwd(i));
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the outgoing edge with label c if it exists and npos otherwise.
 */
edge_index DBG_succ::outgoing_edge_idx(node_index i, TAlphabet c) const {
    CHECK_NODE(i);
    assert(c < alph_size);

    return pick_edge(select_last(i), i, c);
}

/**
 * Given an edge index i and a character c, get the index of the edge with
 * label c outgoing from the same source node if such exists and npos otherwise.
 */
edge_index DBG_succ::pick_edge(edge_index edge, node_index node, TAlphabet c) const {
    CHECK_INDEX(edge);
    CHECK_NODE(node);
    assert(c < alph_size);
    assert(get_source_node(edge) == node);

    uint64_t j = pred_W(edge, c);
    if (!j || get_source_node(j) == node)
        return j;

    j = pred_W(edge, c + alph_size);
    if (!j || get_source_node(j) == node)
        return j;

    return npos;
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the node the edge is pointing to.
 */
node_index DBG_succ::outgoing(node_index i, TAlphabet c) const {
    CHECK_NODE(i);

    c %= alph_size;

    edge_index j = outgoing_edge_idx(i, c);
    if (j == npos)
        return npos;

    uint64_t offset = F_[c];
    uint64_t rank = rank_W(j, c);

    return rank_last(offset) + rank;
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the node the incoming edge belongs to.
 */
node_index DBG_succ::incoming(node_index i, TAlphabet c) const {
    CHECK_NODE(i);

    // only one incoming edge for the dummy source node
    if (i == 1) {
        if (c == kSentinelCode) {
            return 1;
        } else {
            return npos;
        }
    }

    c %= alph_size;

    // check if the first incoming edge has label `c`
    edge_index edge = select_last(i);
    uint64_t x = bwd(edge);

    if (get_minus_k_value(x, k_ - 1).first == c)
        return x ? get_source_node(x) : npos;

    if (x + 1 == W_->size())
        return npos;

    TAlphabet d = get_node_last_value(edge);
    uint64_t y = succ_W(x + 1, d);

    // iterate over the rest of the incoming edges
    while (x + 1 < y) {
        x = succ_W(x + 1, d + alph_size);
        if (x < y && get_minus_k_value(x, k_ - 1).first == c) {
            return x ? get_source_node(x) : npos;
        }
    }
    return npos;
}

DBG_succ::node_index DBG_succ::traverse(node_index node, char edge_label) const {
    return outgoing(node, encode(edge_label));
}

DBG_succ::node_index DBG_succ::traverse_back(node_index node, char edge_label) const {
    return incoming(node, encode(edge_label));
}

void DBG_succ::call_adjacent_incoming_edges(edge_index edge,
                                            std::function<void(edge_index)> callback) const {
    CHECK_INDEX(edge);

    edge_index next_incoming = bwd(edge);

    callback(next_incoming);

    if (next_incoming + 1 == W_->size()) {
        assert(indegree(get_source_node(edge)) == 1);
        return;
    }

    TAlphabet d = get_node_last_value(edge);
    assert(d == get_W(next_incoming));

    // iterate through all indices with edge label d + alph_size
    // which are less than the next index with edge label d
    const auto ubound = succ_W(next_incoming + 1, d);

    while (++next_incoming < ubound
            && (next_incoming = succ_W(next_incoming, d + alph_size)) < ubound) {
        callback(next_incoming);
    }
}

/**
 * Given a node index i, this function returns the number of outgoing
 * edges from node i.
 */
size_t DBG_succ::outdegree(node_index i) const {
    CHECK_NODE(i);

    return select_last(i) - (i == 1 ? 0 : select_last(i - 1));
}

/**
 * Given an edge index i, this function returns true if that is
 * the only outgoing edge from its source node.
 */
bool DBG_succ::is_single_outgoing(edge_index i) const {
    CHECK_INDEX(i);

    return get_last(i) && (i == 1 || get_last(i - 1));
}

/**
 * Given an edge index i, this function returns true if that is
 * the only edge incoming to its target node.
 */
bool DBG_succ::is_single_incoming(edge_index i) const {
    CHECK_INDEX(i);

    TAlphabet c = get_W(i);
    if (c >= alph_size)
        return false;

    // start from the next edge
    i++;

    // trying to avoid calls of succ_W
    size_t max_iter = 1000;
    size_t end = std::min(W_->size(), i + max_iter);

    while (i < end) {
        if (get_W(i) == c + alph_size)
            return false;
        if (get_W(i) == c)
            return true;
        i++;
    }

    return i == W_->size() || succ_W(i, c) <= succ_W(i, c + alph_size);
}

/**
 * Given a node index i, this function returns the number of incoming
 * edges to the node i.
 */
size_t DBG_succ::indegree(node_index i) const {
    CHECK_NODE(i);

    if (i == 1)
        return 1;

    edge_index edge = select_last(i);
    uint64_t x = bwd(edge);
    if (x + 1 == W_->size())
        return 1;

    TAlphabet d = get_node_last_value(edge);

    uint64_t y = succ_W(x + 1, d);
    return 1 + rank_W(y - 1, d + alph_size) - rank_W(x - 1, d + alph_size);
}


/**
 * Given a string str and a maximal number of edit operations
 * max_distance, this function returns all nodes with labels at most
 * max_distance many edits away from str.
 */
std::vector<HitInfo> DBG_succ::index_fuzzy(const std::string &str,
                                           size_t max_distance) const {
    std::vector<HitInfo> result;
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits;
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits2;
    uint64_t rl;
    uint64_t ru;

    // walk through pattern, thereby collecting possible partial matches
    // once the end of the pattern is reached, add match to results

    // init match/mismatch to first pattern position
    TAlphabet s = encode(str[0]);
    assert(alph_size == 6 && "This function is defined for DNA sequences only");
    // TODO: review and test out this function for protein graphs
    for (TAlphabet b = 1; b < 5; ++b) {
        rl = F_[b] + 1 < W_->size()
             ? succ_last(F_[b] + 1)
             : W_->size();
        ru = b + 1 < alph_size
             ? F_[b + 1]
             : W_->size() - 1;
        //std::cout << "pushing: rl " << rl << " ru " << ru << " str_pos 1 max_distance " << (uint64_t) (b != s) << std::endl;
        //std::cout << "s " << s << " b " << b << std::endl;
        std::vector<uint64_t> tmp;
        hits.push({ rl, ru, 1, 1, static_cast<uint64_t>(b != s),
                    std::string(1, decode(b)), tmp });

        // opening/extending a gap in the pattern starting with the first position
        if (max_distance > 0) {
            for (size_t p = 1; p < str.length() - 1; ++p) {
                TAlphabet ss = encode(str[p]);
                if ((p + (b != ss)) > max_distance)
                    break;
                hits.push({ rl, ru, p + 1, 1, p + (b != ss),
                            std::string(p, 'd') + std::string(1, decode(b)), tmp });
                //std::cout << "a) adding '-'" << std::endl;
            }
        }
    }

    // walk through pattern thereby extending all partial hits
    while (hits.size() > 0) {
        while (hits.size() > 0) {
            HitInfo curr_hit(hits.top());
            hits.pop();
            //std::cout << "loaded: rl " << curr_hit.rl << " ru " << curr_hit.ru << " dist " << curr_hit.distance << std::endl;

            if (curr_hit.str_pos < str.length()) {

                // opening/extending a gap in the graph, leaving current pattern position unmatched
                if (curr_hit.distance < max_distance) {
                    hits2.push({ curr_hit.rl, curr_hit.ru, curr_hit.str_pos + 1,
                                 curr_hit.graph_pos, curr_hit.distance + 1,
                                 curr_hit.cigar + 'd', curr_hit.path });
                    //std::cout << "b) " << curr_hit.cigar << " adding '-'" << std::endl;
                }

                s = encode(str[curr_hit.str_pos]);

                // has the number of matches exceeded the node length?
                // there are three possible scenarios for extension of the path:
                //  1) pattern is shorter than the node length --> get an interval of matching nodes
                //  2) pattern length exactly mathces the node length --> there is one correponding node
                //  3) pattern is longer than the node length --> we append to a path
                if (curr_hit.graph_pos >= k_) {
                //    std::cout << "push back tp path " << curr_hit.rl << std::endl;
                    curr_hit.path.push_back(curr_hit.rl);
                }

                // iterate through all possible extensions of current position
                for (TAlphabet b = 1; b < 5; ++b) {
                    if (curr_hit.distance <= max_distance) {

                        // we cannot afford any more mismatches
                        if ((curr_hit.distance + (b != s)) > max_distance)
                            continue;

                        // re-define range of nodes to check for outgoing nodes
                        rl = std::min(succ_W(pred_last(curr_hit.rl - 1) + 1, b),
                                      succ_W(pred_last(curr_hit.rl - 1) + 1, b + alph_size));
                        ru = std::max(pred_W(curr_hit.ru, b),
                                      pred_W(curr_hit.ru, b + alph_size));

                        // the current range in W does not contain our next symbol
                        if (rl >= W_->size() || ru >= W_->size() || rl > ru)
                            continue;

                        // update the SA range with the current symbol b
                        rl = fwd(rl);
                        ru = fwd(ru);

                        // range is empty
                        if ((rl == 0) && (ru == 0))
                            continue;

                        // add hit for extension in next step
                        hits2.push({ rl, ru, curr_hit.str_pos + 1,
                                     curr_hit.graph_pos + 1, curr_hit.distance + (b != s),
                                     curr_hit.cigar + decode(b), curr_hit.path });

                        // opening/extending a gap in the pattern, leaving current graph position unmatched
                        // --> choose any available mismatching next edge
                        if (b != s) {
                            hits2.push({ rl, ru, curr_hit.str_pos,
                                         curr_hit.graph_pos + 1, curr_hit.distance + 1,
                                         curr_hit.cigar + 'i', curr_hit.path });
                        }
                    }
                }
            } else {
                // collect results
                //std::make_pair(curr_hit.rl < curr_hit.ru ? curr_hit.ru : curr_hit.rl, curr_hit.cigar));
                result.push_back(curr_hit);
            }
        }
        hits.swap(hits2);
    }

    return result;
}


/**
 * Given a node label kmer, this function returns the index
 * of the corresponding node or the closest predecessor, if no node
 * with the sequence is not found.
 */
node_index DBG_succ::pred_kmer(const std::vector<TAlphabet> &kmer) const {
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
        assert(s < alph_size);

        uint64_t last_target = std::max(pred_W(last_, s),
                                        pred_W(last_, s + alph_size));
        if (last_target > 0) {
            if (rank_last(last_target - 1) < rank_last(last_ - 1))
                shift = 0;
            last_ = fwd(last_target);
            continue;
        }
        assert(s > 0);

        last_target = std::min(succ_W(last_, s),
                               succ_W(last_, s + alph_size));

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
bool DBG_succ::compare_node_suffix(edge_index first, edge_index second) const {
    for (size_t i = 0; i < k_ - 1; ++i) {
        if (get_node_last_value(first) != get_node_last_value(second)) {
            return false;
        }
        first = bwd(first);
        second = bwd(second);
    }
    return true;
}

bool DBG_succ::compare_node_suffix(edge_index first, const TAlphabet *second) const {
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
std::vector<TAlphabet> DBG_succ::get_node_seq(edge_index k_node) const {
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
std::string DBG_succ::get_node_str(edge_index k_node) const {
    CHECK_INDEX(k_node);
    return decode(get_node_seq(k_node));
}

void DBG_succ::map_to_nodes(const std::string &sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    auto seq_encoded = encode(sequence);

    for (size_t i = 0; i + k_ - 1 < seq_encoded.size(); ++i) {
        auto node = map_to_node(&seq_encoded[i], &seq_encoded[i + k_]);

        callback(node);

        if (terminate())
            return;

        if (!node)
            continue;

        while (i + k_ < seq_encoded.size()) {
            node = outgoing(node, seq_encoded[i + k_]);
            if (!node)
                break;

            callback(node);

            if (terminate())
                return;

            i++;
        }
    }
}

std::vector<node_index> DBG_succ::map_to_nodes(const std::string &sequence,
                                               size_t kmer_size) const {
    if (kmer_size == 0 || kmer_size > k_)
        kmer_size = k_;

    if (sequence.size() < kmer_size)
        return {};

    auto seq_encoded = encode(sequence);

    std::vector<node_index> indices;

    for (size_t i = 0; i + kmer_size <= seq_encoded.size(); ++i) {
        edge_index edge = index_range(&seq_encoded[i],
                                      &seq_encoded[i + kmer_size]).second;
        node_index node = edge ? get_source_node(edge) : npos;
        indices.push_back(node);

        if (!edge || kmer_size != k_ || !indices.back())
            continue;

        // This boost is valid only if alignment length equals k since
        // otherwise, when alignment length is less than k,
        // the existence of an edge between node ending with preceeding
        // aligned substring and node ending with the next aligned substring
        // is not guaranteed. It may be a suffix of some other k-mer.
        while (i + kmer_size < seq_encoded.size()) {
            node = outgoing(node, seq_encoded[i + kmer_size]);
            if (!node)
                break;

            indices.push_back(node);
            i++;
        }
    }

    return indices;
}

void DBG_succ::map_to_edges(const std::string &sequence,
                            const std::function<void(edge_index)> &callback,
                            const std::function<bool()> &terminate) const {
    auto seq_encoded = encode(sequence);

    for (size_t i = 0; i + k_ + 1 <= seq_encoded.size(); ++i) {
        auto edge = map_to_edge(&seq_encoded[i], &seq_encoded[i + k_ + 1]);

        callback(edge);

        if (terminate())
            return;

        while (edge && i + k_ + 1 < seq_encoded.size()) {
            edge = fwd(edge);
            edge = pick_edge(edge, get_source_node(edge), seq_encoded[i + k_ + 1]);

            callback(edge);

            if (terminate())
                return;

            i++;
        }
    }
}

std::vector<edge_index> DBG_succ::map_to_edges(const std::string &sequence) const {
    std::vector<edge_index> indices;
    indices.reserve(sequence.size());
    map_to_edges(sequence,
                 [&indices](edge_index i) { indices.push_back(i); });
    return indices;
}

bool DBG_succ::find(const std::string &sequence,
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

bool DBG_succ::find(const std::string &sequence,
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
        auto edge = map_to_edge(&seq_encoded[i], &seq_encoded[i + kmer_size]);
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
            edge = pick_edge(edge, get_source_node(edge), seq_encoded[i + kmer_size]);
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

        auto edge = map_to_edge(&seq_encoded[i], &seq_encoded[i + kmer_size]);
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
            edge = pick_edge(edge, get_source_node(edge), seq_encoded[i + kmer_size]);
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

std::vector<std::vector<HitInfo>> DBG_succ::align_fuzzy(const std::string &sequence,
                                                        size_t alignment_length,
                                                        size_t max_distance) const {
    std::vector<std::vector<HitInfo>> hit_list;

    if (alignment_length == 0) {

    } else {
        alignment_length = alignment_length < 2 ? 2 : alignment_length;
        alignment_length = alignment_length < sequence.size()
                                    ? alignment_length
                                    : sequence.size();
        for (size_t i = 0; i < sequence.size() - alignment_length + 1; ++i) {
            std::string kmer(sequence.data() + i, sequence.data() + i + alignment_length);
            hit_list.push_back(index_fuzzy(kmer, max_distance));
        }
    }
    return hit_list;
}

/**
 * Returns the number of nodes on the current graph.
 */
uint64_t DBG_succ::num_nodes() const {
    return last_->num_set_bits();
}

/**
 * Return the number of edges in the current graph.
 */
uint64_t DBG_succ::num_edges() const {
    return W_->size() - 1;
}

/**
 * This function gets a value of the alphabet c and updates the offset of
 * all following values by +1 is positive is true and by -1 otherwise.
 */
void DBG_succ::update_F(TAlphabet c, int value) {
    assert(c < alph_size);
    assert(std::abs(value) == 1);

    for (TAlphabet i = c + 1; i < alph_size; i++) {
        F_[i] += value;
    }
}

TAlphabet DBG_succ::encode(char s) const {
    assert(kmer_extractor_.encode(kSentinel) != kSentinelCode);
    assert(kmer_extractor_.encode(kSentinel) != kSentinelCode + alph_size);
    return kmer_extractor_.encode(s);
}

std::vector<TAlphabet> DBG_succ::encode(const std::string &sequence) const {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), [this](char c) { return this->encode(c); });
    return seq_encoded;
}

char DBG_succ::decode(TAlphabet c) const {
    assert(kmer_extractor_.alphabet[kSentinelCode] == kSentinel);
    assert(c < 2 * alph_size);
    return kmer_extractor_.decode(c % alph_size);
}

std::string DBG_succ::decode(const std::vector<TAlphabet> &sequence) const {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(), [this](TAlphabet x) { return this->decode(x); });
    return str;
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

void DBG_succ::switch_state(Config::StateType new_state) {

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
        case Config::DYN: {
            convert<wavelet_tree_dyn, bit_vector_dyn>(&W_, &last_);
            break;
        }
    }
    state = new_state;
}

void DBG_succ::print_internal_representation(std::ostream &os) const {
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

void DBG_succ::print(std::ostream &os) const {
    assert(is_valid());
    auto vertex_header = std::string("Vertex");
    vertex_header.resize(k_, ' ');

    os << "Index" << "\t" << "L"
                  << "\t" << vertex_header
                  << "\t" << "W" << std::endl;

    for (uint64_t i = 1; i < W_->size(); i++) {
        os << i << "\t" << get_last(i)
                << "\t" << get_node_str(i)
                << "\t" << decode(get_W(i))
                        << (get_W(i) >= alph_size
                                ? "-"
                                : "")
                        << std::endl;
    }
}

void DBG_succ::print_adj_list(std::ostream &os) const {
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
void DBG_succ::add_sequence(const std::string &seq, bool try_extend,
                            std::vector<uint64_t> *edges_inserted) {
    if (seq.size() < k_ + 1)
        return;

    auto sequence = encode(seq);

    uint64_t source;

    if (!try_extend || !(source = index(sequence.data(), sequence.data() + k_))) {
        sequence.insert(sequence.begin(), k_, kSentinelCode);
        source = 1; // the dummy source node
    }

    for (size_t i = 0; i < sequence.size() - k_; ++i) {
        // print the process
        if (i > 0 && i % 1'000 == 0) {
            verbose_cout(".");
            if (i % 10'000 == 0)
                verbose_cout(i, " - edges ", num_edges(), " / nodes ", num_nodes(), "\n");
        }

        source = append_pos(sequence[i + k_], source, &sequence[i], edges_inserted);
    }

    verbose_cout("edges ", num_edges(), " / nodes ", num_nodes(), "\n");
}

/**
 * Given a character c and an edge index, this function
 * creates an outgoing edge from the same source node with
 * label c if it is not a part of the graph yet.
 */
edge_index DBG_succ::append_pos(TAlphabet c, edge_index source_node,
                                const TAlphabet *source_node_kmer,
                                std::vector<uint64_t> *edges_inserted) {
    CHECK_INDEX(source_node);
    assert(source_node_kmer);
    assert(std::vector<TAlphabet>(source_node_kmer, source_node_kmer + k_)
                                                == get_node_seq(source_node));

    // get range of identical nodes (without W) pos current end position
    uint64_t begin = pred_last(source_node - 1) + 1;
    uint64_t end = succ_last(source_node) + 1;

    // get position of the first occurence of c or c- in W after p
    uint64_t prev_c_pos = std::max(pred_W(end - 1, c),
                                   pred_W(end - 1, c + alph_size));
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
    last_->insertBit(sentinel_pos, true);

    if (edges_inserted)
        edges_inserted->push_back(sentinel_pos);

    assert((*W_)[0] == 0);

    return sentinel_pos;
}


uint64_t DBG_succ::insert_edge(TAlphabet c, uint64_t begin, uint64_t end) {
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
        last_->insertBit(begin, false);
        W_->insert(pos, c);

        assert(pos);
        return pos;
    }
}


// Given an edge list, remove them from the graph.
// TODO: fix the implementation (anchoring the isolated nodes)
void DBG_succ::erase_edges_dyn(const std::set<edge_index> &edges) {
    uint64_t shift = 0;

    for (edge_index edge : edges) {
        assert(edge >= shift);
        uint64_t edge_id = edge - shift;

        uint64_t d = get_W(edge_id);
        if (d < alph_size) {
            //fix W array
            uint64_t next = edge_id + 1;
            uint64_t j = next < W_->size()
                            ? succ_W(next, d)
                            : W_->size();
            for (uint64_t i = next; i < j; ++i) {
                if (get_W(i) == d + alph_size) {
                    W_->set(i, d);
                    break;
                }
            }
        }
        W_->remove(edge_id);
        update_F(get_node_last_value(edge_id), -1);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        if (get_last(edge_id) && (edge >= shift + 1)
                              && !get_last(edge_id - 1)) {
            last_->deleteBit(edge_id - 1);
        } else {
            last_->deleteBit(edge_id);
        }
        shift++;
    }
}

/**
 * Erase exactly all the masked edges from the graph,
 * may invalidate the graph (if leaves nodes with no incoming edges).
 * Returns the number of edges erased.
 */
uint64_t DBG_succ::erase_edges(const std::vector<bool> &edges_to_remove_mask) {
    size_t num_edges_to_remove = std::count(edges_to_remove_mask.begin(),
                                            edges_to_remove_mask.end(), true);
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
    std::vector<bool> first_removed(alph_size, false);

    for (size_t i = 0, new_i = 0; i < edges_to_remove_mask.size(); ++i) {
        TAlphabet c = get_W(i);
        if (edges_to_remove_mask[i]) {
            if (c < alph_size)
                first_removed[c] = true;
        } else {
            if (c >= alph_size && first_removed[c % alph_size]) {
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
void DBG_succ::edge_DFT(edge_index start,
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
template <typename T, typename U>
void traverse_dummy_edges(const DBG_succ &graph,
                          edge_index subtree_root,
                          size_t check_depth,
                          std::vector<T> *redundant_mask,
                          std::vector<T> *traversed_mask,
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
                traversed_mask->at(edge) = 1;

            if (redundant_path.back())
                redundant_mask->at(edge) = 2;

            redundant_path.pop_back();
        },
        [&](edge_index edge) {
            if (redundant_path.size() == check_depth) {
                if (!redundant_mask
                        || (!redundant_mask->at(edge)
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
uint64_t traverse_dummy_edges(const DBG_succ &graph,
                              std::vector<bool> *redundant_mask,
                              std::vector<bool> *traversed_mask,
                              size_t num_threads,
                              bool verbose) {
    assert(!redundant_mask || redundant_mask->size() == graph.get_W().size());
    assert(!traversed_mask || traversed_mask->size() == graph.get_W().size());

    if (traversed_mask)
        traversed_mask->at(1) = true;

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

    utils::ThreadPool pool(num_threads);

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
        for (size_t i = 0; i < edges_threadsafe->size(); ++i) {
            if ((*edges_threadsafe)[i] && traversed_mask)
                traversed_mask->at(i) = true;
            if ((*edges_threadsafe)[i] == 2 && redundant_mask)
                redundant_mask->at(i) = true;
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
 * Return value: edges removed from the initial graph.
 */
std::vector<bool>
DBG_succ::erase_redundant_dummy_edges(std::vector<bool> *source_dummy_edges,
                                      size_t num_threads,
                                      bool verbose) {
    std::vector<bool> redundant_dummy_edges_mask(W_->size(), false);

    if (source_dummy_edges) {
        source_dummy_edges->assign(W_->size(), false);
        source_dummy_edges->at(1) = true;
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

uint64_t DBG_succ::mark_source_dummy_edges(std::vector<bool> *mask,
                                           size_t num_threads,
                                           bool verbose) const {
    assert(!mask || mask->size() == W_->size());

    return traverse_dummy_edges(*this, NULL, mask, num_threads, verbose);
}

uint64_t DBG_succ::mark_sink_dummy_edges(std::vector<bool> *mask) const {
    if (!mask)
        return rank_W(num_edges(), 0) - 1;

    assert(mask->size() == W_->size());

    uint64_t num_dummy_sink_edges = 0;

    // skip the main dummy source
    for (uint64_t i = 2; i < W_->size(); ++i) {
        assert(get_W(i) != alph_size);
        if (!get_W(i)) {
            mask->at(i) = true;
            num_dummy_sink_edges++;
        }
    }

    assert(num_dummy_sink_edges == rank_W(num_edges(), 0) - 1);

    return num_dummy_sink_edges;
}

std::vector<bool> DBG_succ::mark_all_dummy_edges(size_t num_threads) const {
    std::vector<bool> edge_mask(num_edges() + 1, 0);

    mark_source_dummy_edges(&edge_mask, num_threads);
    mark_sink_dummy_edges(&edge_mask);

    // exclude 0 as the dummy index that denotes not existing k-mers
    edge_mask[0] = true;

    return edge_mask;
}

std::vector<bool> DBG_succ::prune_and_mark_all_dummy_edges(size_t num_threads) {
    std::vector<bool> edge_mask(num_edges() + 1, 0);

    erase_redundant_dummy_edges(&edge_mask, num_threads);
    mark_sink_dummy_edges(&edge_mask);

    // exclude 0 as the dummy index that denotes not existing k-mers
    edge_mask[0] = true;

    return edge_mask;
}

/**
 * Get the nodes of graph |other| merged into the current graph. The graph
 * |other| is fully traversed and all edges are added to to the current graph.
 * This function is well suited to merge small graphs into large ones.
 */
void DBG_succ::merge(const DBG_succ &other) {
    other.call_sequences([&](const std::string &sequence) {
        add_sequence(sequence, true);
    });
}

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void DBG_succ::call_paths(Call<const std::vector<edge_index>,
                               const std::vector<TAlphabet>&> callback,
                          bool split_to_contigs) const {
    // keep track of reached edges
    std::vector<bool> discovered(W_->size(), false);
    // keep track of edges that are already included in covering paths
    std::vector<bool> visited(W_->size(), false);

    // start at all nodes with more than one outgoing edges
    for (uint64_t i = 1; i < W_->size(); ++i) {
        if (!visited[i] && !is_single_outgoing(i))
            call_paths(i, callback, split_to_contigs, &discovered, &visited);
    }

    // process all the cycles left that have not beed traversed
    for (uint64_t i = 1; i < W_->size(); ++i) {
        if (!visited[i])
            call_paths(i, callback, split_to_contigs, &discovered, &visited);
    }
}

struct Edge {
    DBG_succ::edge_index id;
    std::vector<TAlphabet> source_kmer;
};

void DBG_succ::call_paths(edge_index starting_kmer,
                          Call<const std::vector<edge_index>,
                               const std::vector<TAlphabet>&> callback,
                          bool split_to_contigs,
                          std::vector<bool> *discovered_ptr,
                          std::vector<bool> *visited_ptr) const {
    assert(discovered_ptr && visited_ptr);

    auto &discovered = *discovered_ptr;
    auto &visited = *visited_ptr;
    // store all branch nodes on the way
    std::vector<uint64_t> path;
    std::vector<TAlphabet> kmer;

    discovered[starting_kmer] = true;
    std::deque<Edge> edges { { starting_kmer, get_node_seq(starting_kmer) } };

    // keep traversing until we have worked off all branches from the queue
    while (!edges.empty()) {
        uint64_t edge = edges.front().id;
        auto sequence = std::move(edges.front().source_kmer);
        path.clear();
        edges.pop_front();

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (!visited[edge]) {
            assert(edge > 0 && discovered[edge]);

            // visit the edge
            sequence.push_back(get_W(edge) % alph_size);
            path.push_back(edge);
            visited[edge] = true;

            // stop traversing if the next node is a dummy sink
            if (!sequence.back())
                break;

            // stop traversing if we call contigs and this
            // is not the only incoming edge
            bool continue_traversal = !split_to_contigs || is_single_incoming(edge);

            // make one traversal step
            edge = fwd(edge);

            // traverse if there is only one outgoing edge
            if (continue_traversal && is_single_outgoing(edge)) {
                discovered[edge] = true;
                continue;
            }

            kmer.assign(sequence.end() - k_, sequence.end());
            edge_index next_edge = 0;

            // loop over outgoing edges
            do {
                if (!next_edge && !split_to_contigs && !visited[edge]) {
                    // save the edge for visiting if we extract arbitrary paths
                    discovered[edge] = true;
                    next_edge = edge;
                } else if (!discovered[edge]) {
                    // discover other edges
                    discovered[edge] = true;
                    edges.push_back({ edge, kmer });
                }
            } while (--edge > 0 && !get_last(edge));

            // stop traversing this sequence if the next edge was not selected
            if (!next_edge)
                break;

            // pick the last outgoing but not yet visited
            // edge and continue traversing the graph
            edge = next_edge;
        }

        if (path.size())
            callback(path, sequence);
    }
}

void DBG_succ::call_sequences(Call<const std::string&> callback) const {
    std::string sequence;

    call_paths([&](const auto&, const auto &path) {
        sequence.clear();

        for (TAlphabet c : path) {
            if (c != DBG_succ::kSentinelCode) {
                sequence.push_back(DBG_succ::decode(c));
            }
        }

        if (sequence.size())
            callback(sequence);
    }, false);
}

void DBG_succ::call_contigs(Call<const std::string&> callback,
                            size_t max_pruned_dead_end_size) const {
    std::string sequence;

    call_paths([&](const auto&, const auto &path) {
        sequence.clear();

        assert(path.size());

        for (TAlphabet c : path) {
            if (c != kSentinelCode) {
                sequence.push_back(DBG_succ::decode(c));
            }
        }

        if (!sequence.size())
            return;

        if ((path.front() != kSentinelCode && path.back() != kSentinelCode)
                || sequence.size() > k_ + max_pruned_dead_end_size)
            callback(sequence);

    }, true);
}

void DBG_succ::call_edges(Call<edge_index,
                               const std::vector<TAlphabet>&> callback) const {
    call_paths([&](const auto &indices, const auto &path) {
        assert(path.size() == indices.size() + k_);

        for (size_t i = 0; i < indices.size(); ++i) {
            callback(indices[i],
                     std::vector<TAlphabet>(path.begin() + i,
                                            path.begin() + i + k_ + 1));
        }
    });
}

struct Node {
    DBG_succ::node_index id;
    std::string kmer;
};

/**
 * Traverse graph and iterate over all nodes
 */
void DBG_succ::call_kmers(Call<node_index, const std::string&> callback) const {
    // std::vector<bool> discovered(W_->size(), false);
    std::vector<bool> visited(W_->size(), false);

    // store all branch nodes on the way
    std::queue<Node> branchnodes;

    // start at the source node
    for (uint64_t i = 1; i < W_->size(); ++i) {
        if (!get_last(i) || visited[i])
            continue;

        //TODO: traverse backwards

        branchnodes.push({ i, get_node_str(i) });

        // keep traversing until we have worked off all branches from the queue
        while (!branchnodes.empty()) {
            uint64_t node = branchnodes.front().id;
            std::string kmer = std::move(branchnodes.front().kmer);
            branchnodes.pop();

            // traverse forwards until we reach a sink or
            // the first node that has been already visited
            while (!visited[node]) {
                assert(node > 0);

                if (kmer.front() != DBG_succ::kSentinel)
                    callback(node, kmer);

                visited[node] = true;

                // stop traversing if it's a sink
                if (!get_W(node))
                    break;

                std::copy(kmer.begin() + 1, kmer.end(), kmer.begin());

                // traverse if there is only one outgoing edge
                if (is_single_outgoing(node)) {
                    node = fwd(node);
                    kmer.back() = DBG_succ::decode(get_node_last_value(node));
                } else {
                    // loop over outgoing edges
                    do {
                        assert(get_W(node));

                        uint64_t target_node = fwd(node);

                        if (target_node && !visited[target_node]) {
                            kmer.back() = DBG_succ::decode(get_node_last_value(target_node));
                            branchnodes.push({ target_node, kmer });
                        }
                    } while (--node > 1 && !get_last(node));
                    break;
                }
            }
        }
    }
}

bool DBG_succ::is_valid() const {
    assert((*W_)[0] == 0);
    assert(W_->size() >= 2);
    assert(get_node_str(1) == std::string(k_, kSentinel) && "First kmer must be dummy");
    assert(get_W(1) == kSentinelCode && "First kmer must be dummy");

    for (uint64_t i = 1; i < W_->size(); i++) {
        if (get_node_last_value(i) >= alph_size || get_W(i) >= 2 * alph_size)
            return false;

        auto index_pred = bwd(i);
        if (index_pred < 1
                || index_pred >= W_->size()
                || get_node_last_value(index_pred) >= alph_size
                || get_W(index_pred) >= 2 * alph_size)
            return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream &os, const DBG_succ &graph) {
    graph.print(os);
    return os;
}


DBGSuccinct::DBGSuccinct(size_t k, bool canonical_mode)
      : boss_graph_(std::make_unique<DBG_succ>(k - 1)),
        canonical_mode_(canonical_mode) {}

DBGSuccinct::DBGSuccinct(DBG_succ *boss_graph, bool canonical_mode)
      : boss_graph_(boss_graph),
        canonical_mode_(canonical_mode) {}

size_t DBGSuccinct::get_k() const {
    return boss_graph_->get_k() + 1;
}

// Check whether graph contains fraction of nodes from the sequence
bool DBGSuccinct::find(const std::string &sequence,
                       double discovery_fraction) const {
    return boss_graph_->find(sequence, discovery_fraction);
}

// Traverse the outgoing edge
node_index DBGSuccinct::traverse(node_index node, char next_char) const {
    assert(node);

    // dbg node is a boss edge
    DBG_succ::edge_index edge = boss_graph_->fwd(kmer_to_boss_index(node));
    return boss_to_kmer_index(
        boss_graph_->pick_edge(edge,
                               boss_graph_->get_source_node(edge),
                               boss_graph_->encode(next_char))
    );
}

// Traverse the incoming edge
node_index DBGSuccinct::traverse_back(node_index node, char prev_char) const {
    assert(node);

    // map dbg node, i.e. a boss edge, to a boss node
    auto boss_edge = kmer_to_boss_index(node);

    auto boss_node = boss_graph_->get_source_node(boss_edge);
    auto source_node = boss_graph_->traverse_back(boss_node, prev_char);

    if (!source_node)
        return npos;

    return boss_to_kmer_index(
        boss_graph_->outgoing_edge_idx(source_node,
                                       boss_graph_->get_node_last_value(boss_edge))
    );
}

void DBGSuccinct::call_outgoing_kmers(node_index node,
                                      const OutgoingEdgeCallback &callback) const {
    assert(node);

    auto boss_edge = kmer_to_boss_index(node);

    auto last = boss_graph_->fwd(boss_edge);
    auto first = boss_graph_->pred_last(last - 1) + 1;

    for (auto i = first; i <= last; ++i) {
        assert(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size
                == boss_graph_->get_node_last_value(i));

        auto next = boss_to_kmer_index(i);
        if (next != npos)
            callback(next, boss_graph_->decode(boss_graph_->get_W(i)));
    }
}

void DBGSuccinct::adjacent_outgoing_nodes(node_index node,
                                          std::vector<node_index> *target_nodes) const {
    assert(node);
    assert(target_nodes);

    auto boss_edge = kmer_to_boss_index(node);

    auto last = boss_graph_->fwd(boss_edge);
    auto first = boss_graph_->pred_last(last - 1) + 1;

    for (auto i = first; i <= last; ++i) {
        assert(boss_graph_->get_W(boss_edge) % boss_graph_->alph_size
                == boss_graph_->get_node_last_value(i));

        auto next = boss_to_kmer_index(i);
        if (next != npos)
            target_nodes->emplace_back(next);
    }
}

void DBGSuccinct::adjacent_incoming_nodes(node_index node,
                                          std::vector<node_index> *source_nodes) const {
    assert(node);
    assert(source_nodes);

    auto edge = kmer_to_boss_index(node);

    boss_graph_->call_adjacent_incoming_edges(edge,
        [&](DBG_succ::edge_index incoming_boss_edge) {
            assert(boss_graph_->get_W(incoming_boss_edge) % boss_graph_->alph_size
                    == boss_graph_->get_node_last_value(edge));

            auto prev = boss_to_kmer_index(incoming_boss_edge);
            if (prev != npos)
                source_nodes->emplace_back(prev);
        }
    );
}

// Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
// is passed. If passed, |nodes_inserted| must have length equal
// to the number of nodes in graph.
void DBGSuccinct::add_sequence(const std::string &sequence,
                               bit_vector_dyn *nodes_inserted) {
    add_seq(sequence, nodes_inserted);
    if (canonical_mode_) {
        // insert reverse complement sequence as well,
        // to have all canonical k-mers in graph
        std::string sequence_copy = sequence;
        reverse_complement(sequence_copy.begin(), sequence_copy.end());
        add_seq(sequence_copy, nodes_inserted);
    }
}

void DBGSuccinct::add_seq(const std::string &sequence,
                          bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    if (nodes_inserted) {
        std::vector<uint64_t> inserted_indexes;
        inserted_indexes.reserve(sequence.size());

        boss_graph_->add_sequence(sequence, true, &inserted_indexes);

        for (auto i : inserted_indexes) {
            if (valid_edges_.get())
                valid_edges_->insertBit(i, true);
            nodes_inserted->insertBit(boss_to_kmer_index(i), true);
        }
    } else {
        boss_graph_->add_sequence(sequence, true);
    }

    assert(!valid_edges_.get() || !(*valid_edges_)[0]);
}

std::string DBGSuccinct::get_node_sequence(node_index node) const {
    assert(node);
    assert(node <= num_nodes());

    auto boss_edge = kmer_to_boss_index(node);

    return boss_graph_->get_node_str(boss_edge)
            + boss_graph_->decode(boss_graph_->get_W(boss_edge));
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
void DBGSuccinct::map_to_nodes(const std::string &sequence,
                               const std::function<void(node_index)> &callback,
                               const std::function<bool()> &terminate) const {
    if (sequence.size() < get_k())
        return;

    if (canonical_mode_) {
        auto forward = boss_graph_->map_to_edges(sequence);

        std::string sequence_rev_compl = sequence;
        reverse_complement(sequence_rev_compl.begin(), sequence_rev_compl.end());

        auto rev_compl = boss_graph_->map_to_edges(sequence_rev_compl);

        assert(forward.size() == sequence.size() - get_k() + 1);
        assert(forward.size() == rev_compl.size());

        for (size_t i = 0; i < forward.size() && !terminate(); ++i) {
            if (sequence.substr(i, get_k())
                    < sequence_rev_compl.substr(i, get_k())) {
                callback(boss_to_kmer_index(forward[i]));
            } else {
                callback(boss_to_kmer_index(rev_compl[i]));
            }
        }

    } else {
        boss_graph_->map_to_edges(
            sequence,
            [&](DBG_succ::edge_index i) { callback(boss_to_kmer_index(i)); },
            terminate
        );
    }
}

uint64_t DBGSuccinct::num_nodes() const {
    return valid_edges_.get()
                ? valid_edges_->num_set_bits()
                : boss_graph_->num_edges();
}

bool DBGSuccinct::load(const std::string &filename) {
    {
        std::ifstream instream(remove_suffix(filename, kExtension) + kExtension,
                               std::ios::binary);

        if (!boss_graph_->load(instream))
            return false;

        try {
            canonical_mode_ = load_number(instream);
        } catch (...) {
            canonical_mode_ = false;
        }
    }

    // release the old mask
    valid_edges_.reset();

    std::ifstream instream(remove_suffix(filename, kExtension) + kDummyMaskExtension,
                           std::ios::binary);
    if (!instream.good())
        return true;

    // initialize a new vector
    switch (get_state()) {
        case Config::STAT: {
            valid_edges_.reset(new bit_vector_stat());
            break;
        }
        case Config::DYN: {
            valid_edges_.reset(new bit_vector_dyn());
            break;
        }
        case Config::SMALL: {
            valid_edges_.reset(new bit_vector_small());
            break;
        }
    }

    // load the mask of valid edges (all non-dummy including npos 0)
    if (!valid_edges_->load(instream)) {
        std::cerr << "Error: Can't load dummy edge mask." << std::endl;
        return false;
    }

    if (valid_edges_->size() != boss_graph_->num_edges() + 1 || (*valid_edges_)[0]) {
        std::cerr << "Error: Edge mask is not compatible with graph." << std::endl;
        return false;
    }

    return true;
}

void DBGSuccinct::serialize(const std::string &filename) const {
    {
        const auto out_filename = remove_suffix(filename, kExtension) + kExtension;
        std::ofstream outstream(out_filename, std::ios::binary);
        boss_graph_->serialize(outstream);
        serialize_number(outstream, canonical_mode_);

        if (!outstream.good())
            throw std::ios_base::failure("Can't write to file " + out_filename);
    }

    if (!valid_edges_.get())
        return;

    assert((boss_graph_->get_state() == Config::StateType::STAT
                && dynamic_cast<const bit_vector_stat*>(valid_edges_.get()))
        || (boss_graph_->get_state() == Config::StateType::DYN
                && dynamic_cast<const bit_vector_dyn*>(valid_edges_.get()))
        || (boss_graph_->get_state() == Config::StateType::SMALL
                && dynamic_cast<const bit_vector_small*>(valid_edges_.get())));

    const auto out_filename = remove_suffix(filename, kExtension) + kDummyMaskExtension;
    std::ofstream outstream(out_filename, std::ios::binary);
    if (!outstream.good())
        throw std::ios_base::failure("Can't write to file " + out_filename);

    valid_edges_->serialize(outstream);
}

void DBGSuccinct::switch_state(Config::StateType new_state) {
    if (get_state() == new_state)
        return;

    if (valid_edges_.get()) {
        switch (new_state) {
            case Config::STAT: {
                valid_edges_ = std::make_unique<bit_vector_stat>(
                    valid_edges_->convert_to<bit_vector_stat>()
                );
                break;
            }
            case Config::DYN: {
                valid_edges_ = std::make_unique<bit_vector_dyn>(
                    valid_edges_->convert_to<bit_vector_dyn>()
                );
                break;
            }
            case Config::SMALL: {
                valid_edges_ = std::make_unique<bit_vector_small>(
                    valid_edges_->convert_to<bit_vector_small>()
                );
                break;
            }
        }
    }

    boss_graph_->switch_state(new_state);
}

Config::StateType DBGSuccinct::get_state() const {
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::STAT
                || dynamic_cast<const bit_vector_stat*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::DYN
                || dynamic_cast<const bit_vector_dyn*>(valid_edges_.get()));
    assert(!valid_edges_.get()
                || boss_graph_->get_state() != Config::StateType::SMALL
                || dynamic_cast<const bit_vector_small*>(valid_edges_.get()));

    return boss_graph_->get_state();
}

void DBGSuccinct::mask_dummy_kmers(size_t num_threads, bool with_pruning) {
    valid_edges_.reset();

    std::vector<bool> vector = with_pruning
        ? boss_graph_->prune_and_mark_all_dummy_edges(num_threads)
        : boss_graph_->mark_all_dummy_edges(num_threads);

    for (uint64_t i = 0; i < vector.size(); ++i) {
        vector[i] = !vector[i];
    }

    switch (get_state()) {
        case Config::STAT: {
            valid_edges_ = std::make_unique<bit_vector_stat>(std::move(vector));
            break;
        }
        case Config::DYN: {
            valid_edges_ = std::make_unique<bit_vector_dyn>(std::move(vector));
            break;
        }
        case Config::SMALL: {
            valid_edges_ = std::make_unique<bit_vector_small>(std::move(vector));
            break;
        }
    }

    assert(valid_edges_.get());
    assert(valid_edges_->size() == boss_graph_->num_edges() + 1);
    assert(!(*valid_edges_)[0]);
}

uint64_t DBGSuccinct::kmer_to_boss_index(node_index kmer_index) const {
    assert(kmer_index <= num_nodes());

    if (!valid_edges_.get() || !kmer_index)
        return kmer_index;

    return valid_edges_->select1(kmer_index);
}

DBGSuccinct::node_index DBGSuccinct::boss_to_kmer_index(uint64_t boss_index) const {
    assert(boss_index <= boss_graph_->num_edges());
    assert(!valid_edges_.get() || boss_index < valid_edges_->size());

    if (!valid_edges_.get() || !boss_index)
        return boss_index;

    if (!(*valid_edges_)[boss_index])
        return npos;

    return valid_edges_->rank1(boss_index);
}
