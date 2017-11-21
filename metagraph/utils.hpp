#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>
#include <cstdint>
#include <vector>
#include <deque>

#include "dbg_succinct_libmaus.hpp"

namespace utils {

    typedef uint64_t TAlphabet;

    uint64_t kFromFile(std::string infbase);

    /**
    * This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
    * as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
    * returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
    * second value set to true if G2(k2_node) < G1(k1_node).
    */
    std::pair<bool, bool> compare_nodes(DBG_succ *G1, uint64_t k1_node,
                                        DBG_succ *G2, uint64_t k2_node);

    std::pair<std::vector<bool>, uint64_t> compare_nodes(std::vector<DBG_succ*> G,
                                                         std::vector<uint64_t> k,
                                                         std::vector<uint64_t> n,
                                                         size_t &cnt);
    /**
     *  This function checks whether two given strings given as deques are
     *  identical.
     */
    bool compare_seq(std::deque<TAlphabet> s1,
                     std::deque<TAlphabet> s2, size_t start = 0);

    /**
     *  This function checks whether string s1 is lexicographically inverse
     *  greater than s2.
     */
    bool seq_is_greater(std::deque<TAlphabet> s1,
                        std::deque<TAlphabet> s2);

    std::string get_filetype(std::string &fname);
}

#endif // __UTILS_HPP__
