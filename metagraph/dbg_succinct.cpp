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
#include <deque>
#include <algorithm>
#include <string>
#include <cstdio>

#include "dbg_succinct_construct.hpp"

using libmaus2::util::NumberSerialisation;

#define CHECK_INDEX(idx) \
    assert(idx < W->size()); \
    assert(idx > 0)

#if _PROTEIN_GRAPH
const std::string DBG_succ::alphabet = "$ABCDEFGHIJKLMNOPQRSTUVWYZX"
                                       "$ABCDEFGHIJKLMNOPQRSTUVWYZX";
const TAlphabet kCharToNucleotide[128] = {
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
    26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
    16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26,
    26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
    16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26
};
const size_t DBG_succ::kLogSigma = 6;
#elif _DNA_CASE_SENSITIVE_GRAPH
//for case-specific DNA and RNA (U <-> T) data
const std::string DBG_succ::alphabet = "$ACGTNacgt$ACGTNacgt";
const TAlphabet kCharToNucleotide[128] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 6, 5, 7,  5, 5, 5, 8,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  9, 9, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};
const size_t DBG_succ::kLogSigma = 5;
#elif _DNA_GRAPH
//for DNA and RNA (U <-> T) alphabets
const std::string DBG_succ::alphabet = "$ACGTN$ACGTN";
const TAlphabet kCharToNucleotide[128] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};
const size_t DBG_succ::kLogSigma = 4;
#else
static_assert(false,
    "Define an alphabet: either "
    "_DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
);
#endif

static_assert(sizeof(TAlphabet) * 8 >= DBG_succ::kLogSigma,
              "Choose type for TAlphabet properly");

const TAlphabet DBG_succ::alph_size = DBG_succ::alphabet.size() / 2;

const SequenceGraph::node_index SequenceGraph::npos = 0;


DBG_succ::DBG_succ(size_t k)
      : k_(k), last(new bit_vector_dyn()),
               F(alph_size, 0),
               W(new wavelet_tree_dyn(kLogSigma)) {
    last->insertBit(0, false);
    W->insert(0, 0);

    // add the dummy source node
    last->insertBit(1, true);
    W->insert(0, 0);
    for (size_t j = 1; j < alph_size; j++) {
        F[j] = 1;
    }
    assert(is_valid());
}

DBG_succ::DBG_succ(DBGSuccConstructor *builder) : DBG_succ::DBG_succ() {
    assert(builder);

    k_ = builder->get_k();
    builder->build_graph(this);
    assert(is_valid());
}

DBG_succ::~DBG_succ() {
    delete W;
    delete last;
}

/**
 * Given a pointer to a graph structures G1 and G2, the function compares their elements to the
 * each other. It will perform an element wise comparison of the arrays W, last and
 * F and will only check for identity. If any element differs, the function will return
 * false and true otherwise.
 */
bool DBG_succ::equals_internally(const DBG_succ &other) const {
    // compare size
    if (W->size() != other.W->size()) {
        verbose_cout("sizes of graphs differ", "\n",
                     "1: ", W->size(), "\n",
                     "2: ", other.W->size(), "\n");
        return false;
    }

    assert(F.size() == other.F.size());

    // compare last
    for (size_t i = 0; i < W->size(); ++i) {
        if (get_last(i) != other.get_last(i)) {
            verbose_cout("last differs at position ", i, "\n",
                         "1: last[", i, "] = ", get_last(i) , "\n",
                         "2: last[", i, "] = ", other.get_last(i), "\n");
            return false;
        }
    }

    // compare W
    for (size_t i = 0; i < W->size(); ++i) {
        if (get_W(i) != other.get_W(i)) {
            verbose_cout("W differs at position ", i, "\n",
                         "1: W[", i, "] = ", get_W(i) , "\n",
                         "2: W[", i, "] = ", other.get_W(i), "\n");
            return false;
        }
    }

    // compare F
    for (size_t i = 0; i < F.size(); ++i) {
        if (get_F(i) != other.get_F(i)) {
            verbose_cout("F differs at position ", i, "\n",
                         "1: F[", i, "] = ", get_F(i), "\n",
                         "2: F[", i, "] = ", other.get_F(i), "\n");
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

    while (i < W->size() && j < other.W->size()) {

        std::string first_node, second_node;
        bool first_last, second_last;
        char first_label, second_label;

        do {
            first_node = get_node_str(i);
            first_last = get_last(i);
            first_label = DBG_succ::decode(get_W(i));
            i++;
        } while (first_node.find(kSentinel) != std::string::npos && i < W->size());

        do {
            second_node = other.get_node_str(j);
            second_last = other.get_last(j);
            second_label = DBG_succ::decode(other.get_W(j));
            j++;
        } while (second_node.find(kSentinel) != std::string::npos && j < other.W->size());

        if (i == W->size() || j == other.W->size())
            break;

        if (first_node != second_node
                || first_last != second_last
                || first_label != second_label)
            return false;
    }
    return i == W->size() && j == other.W->size();
}

void DBG_succ::serialize(const std::string &outbase) const {

    std::ofstream outstream(outbase + ".dbg");

    // write F values, k, and state
    NumberSerialisation::serialiseNumberVector(outstream, F);
    NumberSerialisation::serialiseNumber(outstream, k_);
    NumberSerialisation::serialiseNumber(outstream, state);
    outstream.flush();

    // write Wavelet Tree
    W->serialise(outstream);
    outstream.flush();

    // write last array
    last->serialise(outstream);
    outstream.close();
}

bool DBG_succ::load(const std::string &infbase) {
    // if not specified in the file, the default for loading is dynamic
    state = Config::DYN;

    try {
        std::ifstream instream(infbase + ".dbg");

        // load F, k, and state
        F = NumberSerialisation::deserialiseNumberVector<uint64_t>(instream);
        k_ = NumberSerialisation::deserialiseNumber(instream);
        state = static_cast<Config::StateType>(NumberSerialisation::deserialiseNumber(instream));

        if (F.size() != alph_size)
            return false;

        // load W and last arrays
        delete W;
        delete last;
        switch (state) {
            case Config::DYN:
                W = new wavelet_tree_dyn(kLogSigma);
                last = new bit_vector_dyn();
                break;
            case Config::STAT:
                W = new wavelet_tree_stat(kLogSigma);
                last = new bit_vector_stat();
                break;
            case Config::SMALL:
                W = new wavelet_tree_small(kLogSigma);
                last = new bit_vector_small();
                break;
        }
        return W->deserialise(instream) && last->deserialise(instream);
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load DBG_succ from "
                  << infbase << "." << std::endl;
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
    assert(i < W->size());

    return i == 0 ? 0 : W->rank(c, i) - (c == 0);
}

/**
 * Uses the array W and gets a count i and a character c from
 * the alphabet and returns the position of the i-th occurence of c in W.
 */
uint64_t DBG_succ::select_W(uint64_t i, TAlphabet c) const {
    assert(i + (c == 0) <= W->rank(c, W->size() - 1));

    return i == 0 ? 0 : W->select(c, i + (c == 0));
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the last index of a character c preceding in W[1..i].
 */
uint64_t DBG_succ::pred_W(uint64_t i, TAlphabet c) const {
    assert(i < W->size());

    if (get_W(i) == c)
        return i;

    return select_W(rank_W(i, c), c);
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the first index of a character c in W[i..N].
 */
uint64_t DBG_succ::succ_W(uint64_t i, TAlphabet c) const {
    assert(i < W->size());

    if (get_W(i) == c)
        return i;

    uint64_t rk = rank_W(i, c);
    if (rk == W->rank(c, W->size() - 1))
        return W->size();

    return select_W(rk + 1, c);
}

/**
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
uint64_t DBG_succ::rank_last(uint64_t i) const {
    assert(i < last->size());

    return i == 0 ? 0 : last->rank1(i);
}

/**
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
uint64_t DBG_succ::select_last(uint64_t i) const {
    assert(i <= last->num_set_bits());

    return i == 0 ? 0 : last->select1(i);
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
uint64_t DBG_succ::pred_last(uint64_t i) const {
    assert(i < last->size());

    return select_last(rank_last(i));
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
uint64_t DBG_succ::succ_last(uint64_t i) const {
    CHECK_INDEX(i);

    if (get_last(i))
        return i;

    uint64_t next_rank = rank_last(i - 1) + 1;

    if (next_rank > last->num_set_bits())
        return last->size();

    return select_last(next_rank);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character.
 */
uint64_t DBG_succ::bwd(uint64_t i) const {
    CHECK_INDEX(i);

    uint64_t node_rank = rank_last(i);
    if (!get_last(i))
        node_rank++;

    // get value of last position in node i
    TAlphabet c = get_node_last_value(i);
    // get the offset for the last position in node i
    uint64_t o = F[c];
    // compute the offset for this position in W and select it
    return select_W(node_rank - rank_last(o), c);
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
    uint64_t o = F[c];
    // get the rank of c in W at position i
    uint64_t r = rank_W(i, c);
    // select the index of the position in last that is rank many positions after offset
    return select_last(rank_last(o) + r);
}


/**
 * Using the offset structure F this function returns the value of the last
 * position of node i.
 */
TAlphabet DBG_succ::get_node_last_value(uint64_t i) const {
    CHECK_INDEX(i);

    if (i == 0)
        return 0;

    for (TAlphabet c = 0; c < alph_size; c++) {
        if (F[c] >= i)
            return c - 1;
    }
    return alph_size - 1;
}

/**
 * Given index i of a node and a value k, this function
 * will return the k-th last character of node i.
 */
std::pair<TAlphabet, uint64_t> DBG_succ::get_minus_k_value(uint64_t i, uint64_t k) const {
    CHECK_INDEX(i);

    for (; k > 0; --k) {
        i = bwd(i);
    }
    return std::make_pair(get_node_last_value(i), bwd(i));
}

/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the outgoing edge with label c.
 */
uint64_t DBG_succ::outgoing_edge_idx(uint64_t i, TAlphabet c) const {
    CHECK_INDEX(i);
    assert(c < alph_size);

    if (i == 0 || i > W->size())
        return 0;

    uint64_t first_pos = pred_last(i - 1) + 1;
    uint64_t last_pos = succ_last(i);

    uint64_t j = std::max(pred_W(last_pos, c),
                          pred_W(last_pos, c + alph_size));

    if (j < first_pos || j >= W->size())
        return 0;

    return j;
}

/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the node the edge is pointing to.
 */
uint64_t DBG_succ::outgoing(uint64_t i, TAlphabet c) const {
    CHECK_INDEX(i);

    c %= alph_size;

    uint64_t j = outgoing_edge_idx(i, c);
    if (j == 0)
        return 0;

    j = fwd(j);
    if (j == 0 || j == W->size())
        return 0;

    return j;
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the node the incoming edge belongs to.
 */
uint64_t DBG_succ::incoming(uint64_t i, TAlphabet c) const {
    CHECK_INDEX(i);

    // no incoming edges for the dummy source
    if (i == 1)
        return 0;

    c %= alph_size;

    // check if the first incoming edge has label `c`
    uint64_t x = bwd(i);

    if (get_minus_k_value(x, k_ - 1).first == c)
        return succ_last(x);

    if (x + 1 == get_W().size())
        return 0;

    TAlphabet d = get_node_last_value(i);
    uint64_t y = succ_W(x + 1, d);

    // iterate over the rest of the incoming edges
    while (x + 1 < y) {
        x = succ_W(x + 1, d + alph_size);
        if (x < y && get_minus_k_value(x, k_ - 1).first == c) {
            return succ_last(x);
        }
    }
    return 0;
}

DBG_succ::node_index DBG_succ::traverse(node_index node, char edge_label) const {
    return outgoing(node, encode(edge_label));
}

DBG_succ::node_index DBG_succ::traverse_back(node_index node, char edge_label) const {
    return incoming(node, encode(edge_label));
}

/**
 * Given a node index i, this function returns the number of outgoing
 * edges from node i.
 */
uint64_t DBG_succ::outdegree(uint64_t i) const {
    CHECK_INDEX(i);

    return (i < W->size()) ? succ_last(i) - pred_last(i - 1) : 0;
}

bool DBG_succ::is_single_outgoing(uint64_t i) const {
    CHECK_INDEX(i);

    return get_last(i) && (i == 1 || get_last(i - 1));
}

/**
 * Given a node index i, this function returns the number of incoming
 * edges to node i.
 */
uint64_t DBG_succ::indegree(uint64_t i) const {
    CHECK_INDEX(i);

    if (i < 2)
        return 0;

    uint64_t x = bwd(i);
    if (x + 1 == W->size())
        return 1;

    TAlphabet d = get_node_last_value(i);

    uint64_t y = succ_W(x + 1, d);
    return 1 + rank_W(y - 1, d + alph_size) - rank_W(x - 1, d + alph_size);
}



/**
 * Given a string str and a maximal number of edit operations
 * max_distance, this function returns all nodes with labels at most
 * max_distance many edits away from str.
 */
std::vector<HitInfo> DBG_succ::index_fuzzy(const std::string &str,
                                           uint64_t max_distance) const {

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
        rl = F[b] + 1 < W->size()
             ? succ_last(F[b] + 1)
             : W->size();
        ru = b + 1 < alph_size
             ? F[b + 1]
             : W->size() - 1;
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
                        if ((rl >= W->size()) || (ru >= W->size()) || (rl > ru))
                            continue;

                        // update the SA range with the current symbol b
                        rl = outgoing(rl, b);
                        ru = outgoing(ru, b);

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
 * Given a node label s, this function returns the index
 * of the corresponding node or the closest predecessor, if no node
 * with the sequence is not found.
 */
uint64_t DBG_succ::pred_kmer(const std::deque<TAlphabet> &kmer) const {
    assert(kmer.size() == k_);

    // get first
    auto kmer_it = kmer.begin();

    uint64_t last = *kmer_it + 1 < alph_size
                    ? F.at(*kmer_it + 1)
                    : W->size() - 1;
    uint64_t shift = 0;

    // update range iteratively while scanning through s
    while (++kmer_it != kmer.end()) {
        TAlphabet s = *kmer_it;
        assert(s < alph_size);

        uint64_t last_target = std::max(pred_W(last, s),
                                        pred_W(last, s + alph_size));
        if (last_target > 0) {
            if (rank_last(last_target - 1) < rank_last(last - 1))
                shift = 0;
            last = succ_last(outgoing(last_target, s));
            continue;
        }
        assert(s > 0);

        last_target = std::min(succ_W(last, s),
                               succ_W(last, s + alph_size));

        if (last_target < W->size()) {
            last = succ_last(outgoing(last_target, s));
            shift = 1;
        } else {
            last = F[s];
            shift = 0;
        }
    }

    assert(pred_last(last - shift) > 0);
    return pred_last(last - shift);
}


/**
 * This function gets two node indices and returns whether the
 * node labels share a k-1 suffix.
 */
bool DBG_succ::compare_node_suffix(uint64_t i1, uint64_t i2) const {
    for (size_t ii = 0; ii < k_ - 1; ++ii) {
        if (get_node_last_value(i1) != get_node_last_value(i2)) {
            return false;
        }
        i1 = bwd(i1);
        i2 = bwd(i2);
    }
    return true;
}

bool DBG_succ::compare_node_suffix(TAlphabet *ref, uint64_t i2) const {
    TAlphabet *i1 = &ref[k_ - 1];
    for (size_t ii=0; ii < k_ - 1; ii++) {
        if (*i1 != get_node_last_value(i2)) {
            return false;
        }
        i1 = &ref[k_ - 2 - ii];
        i2 = bwd(i2);
    }
    return true;
}

/**
 * This function returns true if node i is a terminal node.
 */
bool DBG_succ::is_terminal_node(uint64_t i) const {
    CHECK_INDEX(i);

    return num_nodes() <= 1 || (i > 1 && get_W(i) == kSentinelCode);
}

/**
* Given a node index k_node, this function returns the k-mer sequence of the
* node in a deque data structure.
*/
std::deque<TAlphabet> DBG_succ::get_node_seq(uint64_t k_node) const {
    CHECK_INDEX(k_node);

    std::deque<TAlphabet> ret(k_, get_node_last_value(k_node));

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
std::string DBG_succ::get_node_str(uint64_t k_node) const {
    CHECK_INDEX(k_node);

    std::string node_string(k_, 0);

    auto node_encoding = get_node_seq(k_node);

    std::transform(node_encoding.begin(), node_encoding.end(),
                   node_string.begin(), decode);
    return node_string;
}

void DBG_succ::align(const std::string &sequence,
                     const std::function<void(edge_index)> &callback,
                     const std::function<bool()> &terminate) const {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), encode);

    for (uint64_t i = 0; i < seq_encoded.size() - k_ + 1; ++i) {
        auto range = index_range(&seq_encoded[i],
                                 &seq_encoded[i + k_]);
        auto kmer_index = std::max(range.first, range.second);

        callback(kmer_index);

        if (terminate())
            return;

        if (kmer_index == 0)
            continue;

        while (i < seq_encoded.size() - k_) {
            auto next = outgoing(kmer_index, seq_encoded[i + k_]);
            if (next == 0)
                break;

            kmer_index = next;

            callback(kmer_index);

            if (terminate())
                return;

            i++;
        }
    }
}

bool DBG_succ::find(const std::string &sequence,
                    double kmer_discovery_fraction) const {
    if (sequence.length() < k_)
        return false;

    const size_t num_kmers = sequence.length() - k_ + 1;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    const size_t max_kmers_missing = num_kmers * (1 - kmer_discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;

    align(sequence,
        [&](uint64_t i) {
            if (i > 0) {
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

std::vector<uint64_t> DBG_succ::index(const std::string &sequence,
                                      uint64_t alignment_length) const {
    if (alignment_length == 0 || alignment_length > k_)
        alignment_length = k_;

    if (sequence.size() < alignment_length)
        return {};

    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), encode);

    std::vector<uint64_t> indices;

    for (uint64_t i = 0; i < seq_encoded.size() - alignment_length + 1; ++i) {
        auto range = index_range(seq_encoded.data() + i,
                                 seq_encoded.data() + i + alignment_length);
        indices.push_back(std::max(range.first, range.second));

        if (alignment_length != k_ || !indices.back())
            continue;

        // This boost is valid only if alignment length equals k since
        // otherwise, when alignment length is less than k,
        // the existence of an edge between node ending with preceeding
        // aligned substring and node ending with the next aligned substring
        // is not guaranteed. It may be a suffix of some other k-mer.
        while (i < seq_encoded.size() - alignment_length) {
            auto next = outgoing(indices.back(),
                                 seq_encoded[alignment_length + i]);
            if (next == 0)
                break;

            indices.push_back(next);
            i++;
        }
    }

    return indices;
}

std::vector<std::vector<HitInfo>> DBG_succ::align_fuzzy(const std::string &sequence,
                                                        uint64_t alignment_length,
                                                        uint64_t max_distance) const {
    std::vector<std::vector<HitInfo>> hit_list;

    if (alignment_length == 0) {

    } else {
        alignment_length = alignment_length < 2 ? 2 : alignment_length;
        alignment_length = alignment_length < sequence.size()
                                    ? alignment_length
                                    : sequence.size();
        for (uint64_t i = 0; i < sequence.size() - alignment_length + 1; ++i) {
            std::string kmer(sequence.data() + i, sequence.data() + i + alignment_length);
            hit_list.push_back(this->index_fuzzy(kmer, max_distance));
        }
    }
    return hit_list;
}

/**
 * Returns the number of nodes on the current graph.
 */
uint64_t DBG_succ::num_nodes() const {
    return last->num_set_bits();
}

/**
 * Return the number of edges in the current graph.
 */
uint64_t DBG_succ::num_edges() const {
    return W->size() - 1;
}

/**
 * This function gets a value of the alphabet c and updates the offset of
 * all following values by +1 is positive is true and by -1 otherwise.
 */
void DBG_succ::update_F(TAlphabet c, int value) {
    assert(c < alph_size);
    assert(std::abs(value) == 1);

    for (TAlphabet i = c + 1; i < alph_size; i++) {
        F[i] += value;
    }
}

/**
 * This function gets a local range in W from lower bound l
 * to upper bound u and swaps the inserted element to the
 * righ location.
 */
//TODO: this function can be improved
//TODO: fix_order_in_W_range
void DBG_succ::sort_W_locally(uint64_t l, uint64_t u) {
    assert(l < W->size());
    assert(u < W->size());
    assert(state == Config::DYN);

    for (uint64_t s = u; s > l; --s) {
        auto first = get_W(s - 1);
        auto second = get_W(s);
        if ((second % alph_size) < (first % alph_size)) {
            W->set(s - 1, second);
            W->set(s, first);
        }
    }
    for (uint64_t s = l; s < u; ++s) {
        auto first = get_W(s);
        auto second = get_W(s + 1);
        if ((first % alph_size) > (second % alph_size)) {
            W->set(s + 1, first);
            W->set(s, second);
        }
    }
}

TAlphabet DBG_succ::encode(char s) {
    assert(kCharToNucleotide[kSentinel] != kSentinelCode);
    assert(kCharToNucleotide[kSentinel] != kSentinelCode + alph_size);
    assert(static_cast<size_t>(s) < 128);
    return kCharToNucleotide[static_cast<size_t>(s)];
}

char DBG_succ::decode(TAlphabet c) {
    assert(alphabet[kSentinelCode] == kSentinel);
    assert(alphabet[kSentinelCode + alph_size] == kSentinel);
    assert(c < alphabet.size());
    return alphabet[c];
}

template <class WaveletTree, class BitVector>
void convert(wavelet_tree **W, bit_vector **last) {
    wavelet_tree *W_new = new WaveletTree((*W)->convert_to<WaveletTree>());
    delete *W;
    *W = W_new;

    bit_vector *last_new = new BitVector((*last)->convert_to<BitVector>());
    delete *last;
    *last = last_new;
}

void DBG_succ::switch_state(Config::StateType new_state) {

    //std::cerr << "switching state from " << this->state << " to " << state << std::endl;
    if (state == new_state)
        return;

    switch (new_state) {
        case Config::STAT: {
            convert<wavelet_tree_stat, bit_vector_stat>(&W, &last);
            break;
        }
        case Config::SMALL: {
            convert<wavelet_tree_small, bit_vector_small>(&W, &last);
            break;
        }
        case Config::DYN: {
            convert<wavelet_tree_dyn, bit_vector_dyn>(&W, &last);
            break;
        }
    }
    state = new_state;
}

void DBG_succ::print_state(std::ostream &os) const {
    assert(is_valid());
    auto vertex_header = std::string("Vertex");
    vertex_header.resize(k_, ' ');

    os << "Index" << "\t" << "L"
                  << "\t" << vertex_header
                  << "\t" << "W" << std::endl;

    for (uint64_t i = 1; i < W->size(); i++) {
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
    for (uint64_t edge = 1; edge < W->size(); ++edge) {
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
                            bit_vector_dyn *edges_inserted) {
    assert(!edges_inserted || edges_inserted->size() == W->size());

    if (seq.size() < k_)
        return;

    std::vector<TAlphabet> sequence(seq.size());
    std::transform(seq.begin(), seq.end(), sequence.begin(), encode);

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
 * This function takes a character c and appends it to the end of the graph
 * sequence given that the corresponding node is not a part of the graph yet.
 */
uint64_t DBG_succ::append_pos(TAlphabet c, uint64_t source_node, TAlphabet *ckmer,
                              bit_vector_dyn *edges_inserted) {
    CHECK_INDEX(source_node);

    // get range of identical nodes (without W) pos current end position
    uint64_t begin = pred_last(source_node - 1) + 1;
    uint64_t end = succ_last(source_node) + 1;

    // get position of the first occurence of c or c- in W after p
    uint64_t prev_c_pos = std::max(pred_W(end - 1, c),
                                   pred_W(end - 1, c + alph_size));
    // if the character already exists return its index
    if (prev_c_pos >= begin)
        return fwd(prev_c_pos);

    /**
     * We found that c does not yet exist in the current range and now have to
     * figure out if we need to add c or c- to the range.
     * To do this, we check if there is a preceding position j1 with W[j1] == c
     * whose node shares a k-1 suffix with the current node.
     * If yes, we add c- instead of c.
     */

    bool is_first_incoming = true;
    if (prev_c_pos > 0)
        is_first_incoming = ckmer ? !compare_node_suffix(ckmer, prev_c_pos)
                                  : !compare_node_suffix(begin, prev_c_pos);

    if (!is_first_incoming) {
        // insert the edge
        insert_edge(c + alph_size, begin, end, edges_inserted);
        return fwd(prev_c_pos);
    }

    // adding a new node can influence one of the following nodes sharing the k-1 suffix
    // get position of the first occurence of c after p (including p + 1)
    uint64_t first_c = end < W->size()
                       ? succ_W(end, c)
                       : W->size();

    bool the_only_incoming = true;
    if (first_c < W->size()) {
        the_only_incoming = ckmer != NULL ? !compare_node_suffix(ckmer, first_c)
                                          : !compare_node_suffix(begin, first_c);
        if (!the_only_incoming) {
            // the inserted edge will not be the first incoming for the target node
            // need to adapt the respective cc to a cc-
            W->set(first_c, c + alph_size);
        }
    }

    // insert the edge
    bool shift = insert_edge(c, begin, end, edges_inserted);

    // Add sentinel if the target node is the new dead-end
    if (!the_only_incoming)
        return fwd(first_c + shift);

    uint64_t sentinel_pos = select_last(rank_last(F[c]) + rank_W(begin - 1, c)) + 1;

    update_F(c, +1);
    W->insert(sentinel_pos, kSentinelCode);
    last->insertBit(sentinel_pos, true);
    if (edges_inserted)
        edges_inserted->insertBit(sentinel_pos, true);
    return sentinel_pos;
}


bool DBG_succ::insert_edge(TAlphabet c, uint64_t begin, uint64_t end,
                           bit_vector_dyn *edges_inserted) {
    if (begin > 1 && get_W(begin) == kSentinelCode) {
        // the source node is the dead-end with outgoing sentinel
        // replace this sentinel with proper label
        W->set(begin, c);
        return 0;
    } else {
        // the source node already has some outgoing edges
        update_F(get_node_last_value(begin), +1);
        W->insert(begin, c);
        last->insertBit(begin, false);
        // FYI: this works only for the node based graph colorings
        // TODO: fix this if we annotate edges
        if (edges_inserted)
            edges_inserted->insertBit(begin, true);
        sort_W_locally(begin, end);
        return 1;
    }
}


// Given an edge list, remove them from the graph.
void DBG_succ::remove_edges(const std::set<uint64_t> &edges) {
    uint64_t shift = 0;

    for (const auto &edge : edges) {
        assert(edge >= shift);
        uint64_t edge_id = edge - shift;

        uint64_t d = get_W(edge_id);
        if (d < alph_size) {
            //fix W array
            uint64_t next = edge_id + 1;
            uint64_t j = succ_W(next, d);
            for (uint64_t i = next; i < j; ++i) {
                if (get_W(i) == d + alph_size) {
                    W->set(i, d);
                    break;
                }
            }
        }
        W->remove(edge_id);
        update_F(get_node_last_value(edge_id), -1);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        if (get_last(edge_id) && (edge >= shift + 1)
                              && !get_last(edge_id - 1)) {
            last->deleteBit(edge_id - 1);
        } else {
            last->deleteBit(edge_id);
        }
        shift++;
    }
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

struct Edge {
    DBG_succ::edge_index id;
    std::deque<TAlphabet> source_kmer;
};

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void DBG_succ::call_paths(const PathCallback &callback) const {
    // keep track of reached edges
    std::vector<bool> discovered(W->size(), false);
    // keep track of edges that are already included in covering paths
    std::vector<bool> visited(W->size(), false);
    // store all branch nodes on the way
    std::queue<Edge> edges;
    std::vector<uint64_t> path;

    // start at the source node
    for (uint64_t i = 1; i < W->size(); ++i) {
        if (visited[i])
            continue;

        //TODO: traverse backwards

        discovered[i] = true;
        edges.push({ i, get_node_seq(i) });

        // keep traversing until we have worked off all branches from the queue
        while (!edges.empty()) {
            uint64_t edge = edges.front().id;
            auto sequence = std::move(edges.front().source_kmer);
            path.clear();
            edges.pop();

            // traverse simple path until we reach its tail or
            // the first edge that has been already visited
            while (!visited[edge]) {
                assert(edge > 0 && discovered[edge]);

                sequence.push_back(get_W(edge) % alph_size);
                path.push_back(edge);
                visited[edge] = true;

                // stop traversing if the next node is a dummy sink
                if (!sequence.back())
                    break;

                // make one traversal step
                edge = fwd(edge);

                // traverse if there is only one outgoing edge
                if (is_single_outgoing(edge)) {
                    discovered[edge] = true;
                    continue;
                } else {
                    std::deque<TAlphabet> kmer(sequence.end() - k_, sequence.end());

                    // loop over outgoing edges
                    do {
                        if (!discovered[edge]) {
                            discovered[edge] = true;
                            edges.push({ edge, kmer });
                        }
                    } while (--edge > 0 && !get_last(edge));
                    break;
                }
            }

            if (path.size())
                callback(path, sequence);
        }
    }
}

void DBG_succ::call_sequences(const SequenceCallback &callback) const {
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
    });
}

void DBG_succ::call_edges(const EdgeCallback &callback) const {
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
void DBG_succ::call_kmers(const KmerCallback &callback) const {
    // std::vector<bool> discovered(W->size(), false);
    std::vector<bool> visited(W->size(), false);

    // store all branch nodes on the way
    std::queue<Node> branchnodes;

    // start at the source node
    for (uint64_t i = 1; i < W->size(); ++i) {
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
    assert(W->size() >= 2);
    assert(get_node_str(1) == std::string(k_, kSentinel) && "First kmer must be dummy");
    assert(get_W(1) == kSentinelCode && "First kmer must be dummy");

    for (uint64_t i = 1; i < W->size(); i++) {
        if (get_node_last_value(i) >= alph_size || get_W(i) >= alphabet.size())
            return false;

        auto index_pred = bwd(i);
        if (index_pred < 1
                || index_pred >= W->size()
                || get_node_last_value(index_pred) >= alph_size
                || get_W(index_pred) >= alphabet.size())
            return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream &os, const DBG_succ &graph) {
    graph.print_state(os);
    return os;
}
