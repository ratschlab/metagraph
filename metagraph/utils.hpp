#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>
#include <cstdint>
#include <vector>
#include <deque>

#include "dbg_succinct.hpp"


namespace utils {

    uint64_t kFromFile(const std::string &infbase);

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
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    std::string get_filetype(const std::string &fname);

} // namespace utils

#endif // __UTILS_HPP__
