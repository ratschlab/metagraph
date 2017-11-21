#ifndef __TRAVERSE_HPP__
#define __TRAVERSE_HPP__

#include <stack>
#include <cstdint>
#include <string>
#include <map>
#include <vector>

#include "dbg_succinct_libmaus.hpp"

namespace traverse {

    typedef uint64_t TAlphabet;

    /**
     * This object collects information about branches during graph traversal, so
     * we know where to jump back to when we reached a dead end.
     */
    struct BranchInfo;

    /**
     * This will hold the graph edges that will be written to the SQL graph output.
     */
    struct JoinInfo;


    /**
     * This is a convenience function that pops the last branch and updates the traversal state.
     */
    BranchInfo pop_branch(std::stack<BranchInfo> &branchnodes,
                          uint64_t &seqId, uint64_t &seqPos, uint64_t &nodeId,
                          uint64_t &lastEdge, bool &isFirst);

    bool finish_sequence(DBG_succ* G, std::string &sequence, uint64_t seqId, std::ofstream &SQLstream);

    size_t traverseGraph(DBG_succ* G, std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream);

    void allelesFromSeq(DBG_succ* G, kstring_t &seq, unsigned int f, std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream, bool isRefRun = false, size_t seqNum = 0);

    /**
     * Take the current graph content and return it in SQL
     * format (GA4GH Spec).
     *
     * We will perform one depth first search traversal of the graph. While we will record
     * one long reference string, we will output all sidepaths on the way.
     */
    void toSQL(DBG_succ *G);

}

#endif // __TRAVERSE_HPP__
