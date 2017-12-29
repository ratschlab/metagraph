#include "traverse.hpp"

#include <stack>
#include <cstdint>
#include <string>
#include <map>
#include <vector>

#include <zlib.h>
#include <libmaus2/digest/md5.hpp>
#include <htslib/kseq.h>
KSEQ_INIT(gzFile, gzread)

#include "dbg_succinct.hpp"


namespace traverse {

#ifdef DBGDEBUG
    bool debug = true;
#else
    bool debug = false;
#endif

/**
 * This object collects information about branches during graph traversal, so 
 * we know where to jump back to when we reached a dead end.
 */
struct BranchInfo {
    uint64_t nodeId = 0;
    uint64_t seqId = 0;
    uint64_t seqPos = 0;
    DBG_succ::TAlphabet lastEdge = 0;
};


/**
 * This will hold the graph edges that will be written to the SQL graph output.
 */
struct JoinInfo {
    uint64_t seqId1 = 0;
    uint64_t seqPos1 = 0;
    uint64_t seqId2 = 0;
    uint64_t seqPos2 = 0;
};


/**
 * This object collects information about branches during graph traversal, so
 * we know where to jump back to when we reached a dead end.
 */
struct BranchInfo;

/**
 * This will hold the graph edges that will be written to the SQL graph output.
 */
struct JoinInfo;


size_t traverseGraph(DBG_succ *G,
                     std::vector<JoinInfo> &joins,
                     std::map<std::pair<uint64_t, DBG_succ::TAlphabet>, uint64_t> &branchMap,
                     const std::string &sqlfbase,
                     std::ofstream &SQLstream);

void allelesFromSeq(DBG_succ *G,
                    kstring_t &seq,
                    unsigned int f,
                    std::vector<JoinInfo> &joins,
                    std::map<std::pair<uint64_t, DBG_succ::TAlphabet>, uint64_t> &branchMap,
                    std::ofstream &SQLstream,
                    bool isRefRun = false,
                    size_t seqNum = 0);


/**
 * This is a convenience function that pops the last branch and updates the traversal state.
 */
BranchInfo pop_branch(std::stack<BranchInfo> &branchnodes,
                      uint64_t &seqId, uint64_t &seqPos,
                      uint64_t &nodeId, uint64_t &lastEdge, bool &isFirst) {
    BranchInfo branch = branchnodes.top();
    branchnodes.pop();
    isFirst = true;
    seqPos = branch.seqPos;
    seqId = branch.seqId;
    lastEdge = branch.lastEdge;
    nodeId = branch.nodeId;

    return branch;
}


bool finish_sequence(std::string &sequence,
                     uint64_t seqId,
                     const std::string &sqlfbase,
                     std::ofstream &SQLstream) {
    if (sequence.length() > 0) {
        if (seqId == 1)
            SQLstream << "INSERT INTO FASTA VALUES (1, '"
                                    << sqlfbase << ".fa');" << std::endl;
        std::ofstream stream;
        if (seqId == 1)
            stream.open((sqlfbase + ".fa").c_str());
        else
            stream.open((sqlfbase + ".fa").c_str(), std::ofstream::app);
        stream << ">seq" << seqId << std::endl;
        uint64_t i = 0;
        while ((i + 80) < sequence.length()) {
            stream << sequence.substr(i, 80) << std::endl;
            i += 80;
        }
        if (i != sequence.length())
            stream << sequence.substr(i) << std::endl;
        stream.close();

        if (debug)
            std::cout << sequence << std::endl;

        std::string md5;
        std::ostringstream test;
        test << sequence;
        libmaus2::util::MD5::md5(test.str(), md5);
        SQLstream << "INSERT INTO Sequence VALUES ("
                  << seqId << ", 1, 'seq"
                  << seqId << "', '"
                  << md5 << "', "
                  << sequence.length() << ");" << std::endl;
        sequence.clear();
        return true;
    } else {
        return false;
    }
}


size_t traverseGraph(DBG_succ *G,
                     std::vector<JoinInfo> &joins,
                     std::map<std::pair<uint64_t, DBG_succ::TAlphabet>, uint64_t> &branchMap,
                     const std::string &sqlfbase,
                     std::ofstream &SQLstream) {
    // store all branch nodes on the way
    std::stack<BranchInfo> branchnodes;
    std::map<uint64_t, std::pair<uint64_t, uint64_t> > nodeId2seqId;
    // bool vector that keeps track of visited nodes
    std::vector<bool> visited(G->get_W().size());
    for (std::vector<bool>::iterator it = visited.begin(); it != visited.end(); ++it) {
        *it = false;
    }
    std::string sequence;
    // for nodes with indegree > 1 we store sequence and index of the 
    // sequence that visited them, so we know where to anchor branches into it
    std::map<uint64_t, std::pair<uint64_t, uint64_t> > nodeId2seqPos; 
    // some initializations
    uint64_t nodeId = 1; // start at source node
    // position in currently traversed sequence,
    // will increase with every visited node and be reset upon new branch
    uint64_t seqPos = 0;
    uint64_t seqId = 1;  // first sequence ID is 1
    size_t seqCnt = 1; // number of total sequences
    bool isFirst = true; // is this the first node in a new sequence?
    size_t out = G->outdegree(nodeId);
    std::pair<uint64_t, uint64_t> branchPos;
    BranchInfo branch;
    DBG_succ::TAlphabet val;
    DBG_succ::TAlphabet lastEdge = 0;
    DBG_succ::TAlphabet joinEdge = 0;
    bool joinOpen = false;
    JoinInfo currJoin;

    while (out > 0 || branchnodes.size() > 0) {

        //fprintf(stderr, "1: nodeId %lu out %lu\n", nodeId, out);
        // we have reached the sink but there are unvisited nodes left on the stack
        if (out == 0) {
            if (branchnodes.size() == 0)
                break;
            // get new branch
            branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
            out = G->outdegree(nodeId);
            if (debug)
                fprintf(stderr, " -- popped %" PRIu64 " -- ", nodeId); 
            joinOpen = true;
            joinEdge = lastEdge + 1;
        }

        // we have not visited that node before
        if (!visited.at(nodeId)) {
            visited.at(nodeId) = true;
            seqPos += isFirst ? 0 : 1;
            isFirst = false;
            val = G->get_node_last_char(nodeId);
            sequence.append(1, G->decode(val % G->alph_size));
            // store seq position of this node (we will join to it later)
            if (G->indegree(nodeId) > 1) {
                nodeId2seqPos.insert(std::make_pair(nodeId, std::make_pair(seqId, seqPos)));
            }
        }

        // there is only one child
        if (out == 1) {
            uint64_t next = G->fwd(nodeId);
            // the next node is new
            if (!visited.at(next)) {
                if (joinOpen) {
                    if (sequence.length() > 0) {
                        finish_sequence(sequence, seqCnt++, sqlfbase, SQLstream);
                    }
                    seqId = seqCnt;
                    seqPos = 0;
                    branchPos = nodeId2seqId[nodeId];
                    joins.push_back({ branchPos.first, branchPos.second, seqId, seqPos });
                    joinOpen = false;
                    branchMap.insert(
                        std::make_pair(
                            std::make_pair(nodeId, joinEdge),
                            joins.size() - 1
                        )
                    );
                }
                nodeId = next;
                lastEdge = 0;
            // we have seen the next node before
            } else {
                // look up the sequence info of that node
                joins.push_back({ seqId, seqPos, nodeId2seqPos[next].first,
                                                 nodeId2seqPos[next].second });
                branchMap.insert(
                    std::make_pair(
                        std::make_pair(nodeId, 1ul),
                        joins.size() - 1
                    )
                );
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
                out = G->outdegree(nodeId);
                if (debug)
                    fprintf(stderr, " -- popped %" PRIu64 " -- ", nodeId); 
                joinOpen = true;
                joinEdge = lastEdge + 1;
            }
            if (debug)
                fprintf(stderr, " new nodeId: %" PRIu64 "\n", nodeId);
        // there are several children
        } else {
            size_t cnt = 0;
            bool updated = false;
            for (DBG_succ::TAlphabet c = 1; c < G->alph_size; ++c) {
                uint64_t next = G->outgoing(nodeId, c);
                if (next > 0) {
                    cnt++;
                    // we already handled this edge erlier
                    if (cnt <= lastEdge)
                        continue;

                    lastEdge++;
                    if (!visited.at(next)) {
                        // there are remaining branches - push node to stack
                        if (cnt < out && next != nodeId) {
                            // each node can branch from exactly one sequence ID
                            if (nodeId2seqId.find(nodeId) == nodeId2seqId.end())
                                nodeId2seqId.insert(
                                    std::make_pair(
                                        nodeId,
                                        std::make_pair(seqId, seqPos)
                                    )
                                );
                            // push node information to stack
                            branchnodes.push({ nodeId, seqId, seqPos, lastEdge });
                            if (debug)
                                fprintf(stderr, " -- pushed %" PRIu64
                                                " : seqId %" PRIu64
                                                " seqPos %" PRIu64
                                                " lastEdge %" PRIu64 " -- ",
                                                nodeId, seqId, seqPos, lastEdge); 
                        }
                        if (joinOpen) {
                            if (sequence.length() > 0) {
                                finish_sequence(sequence, seqCnt++, sqlfbase, SQLstream);
                            }
                            seqId = seqCnt;
                            seqPos = 0;
                            branchPos = nodeId2seqId[nodeId];
                            joins.push_back({ branchPos.first,
                                              branchPos.second, seqId, seqPos });
                            joinOpen = false;
                            branchMap.insert(
                                std::make_pair(
                                    std::make_pair(nodeId, lastEdge),
                                    joins.size() - 1
                                )
                            );
                        }

                        nodeId = next;
                        updated = true;
                        lastEdge = 0;
                        break;
                    } else {
                        // look up the sequence info of that node
                        if (nodeId == next) { 
                            joins.push_back({ nodeId2seqPos[next].first, nodeId2seqPos[next].second,
                                              nodeId2seqPos[next].first, nodeId2seqPos[next].second });
                        } else {
                            joins.push_back({ seqId, seqPos, nodeId2seqPos[next].first,
                                                             nodeId2seqPos[next].second});
                        }
                        branchMap.insert(
                            std::make_pair(
                                std::make_pair(nodeId, lastEdge),
                                joins.size() - 1
                            )
                        );
                    }
                }
            }
            // we are done with this branch
            // we should end up here, when nodes branch to themselves with their last edge
            if (!updated) {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
                out = G->outdegree(nodeId);
                if (debug)
                    fprintf(stderr, " -- popped %" PRIu64 " -- ", nodeId); 
                joinOpen = true;
                joinEdge = lastEdge + 1;
            }
            //if (debug)
            //    fprintf(stderr, " new nodeId: %lu\n", nodeId);
        }
        out = G->outdegree(nodeId);
    }
    // for completeness
    if (sequence.length() > 0)
        finish_sequence(sequence, seqCnt++, sqlfbase, SQLstream);
    else
        seqCnt--;

    if (debug) {
        // output joins
        for (size_t i = 0; i < joins.size(); ++i) {
            std::cout << "(" << joins.at(i).seqId1 << ":"
                             << joins.at(i).seqPos1 << "--"
                             << joins.at(i).seqId2 << ":"
                             << joins.at(i).seqPos2 << ")" << std::endl;
        }
        // output branches for tracking
        for (auto it = branchMap.begin(); it != branchMap.end(); ++it) {
            std::cout << "["
                        << it->first.first
                        << ", "
                        << it->first.second
                      << "] -- " << it->second << std::endl;
        }
    }

    return seqCnt;
}


void allelesFromSeq(DBG_succ* G, kstring_t &seq, unsigned int f,
                    std::vector<JoinInfo> &joins,
                    std::map< std::pair<uint64_t, DBG_succ::TAlphabet>, uint64_t > &branchMap,
                    std::ofstream &SQLstream,
                    bool isRefRun, size_t seqNum) {
    
    uint64_t nodeId = 1;
    uint64_t seqId = 1;
    uint64_t seqPos = 0;
    uint64_t alleleSeqPos = 0;
    // iterate over nodes
    uint64_t out = 0;
    uint64_t edge = 0;
    uint64_t currStart = 0;
    unsigned int alleleCnt = 1;
    DBG_succ::TAlphabet seqVal, seqValNext;
    DBG_succ::TAlphabet nodeVal;
    JoinInfo currJoin;
    std::vector<bool> isRef;

    if (!isRefRun) {
        SQLstream << "INSERT INTO Allele VALUES (" << f+1 << ", 1, NULL);" << std::endl;
    } else {
        for (size_t i = 0; i < seqNum; ++i)
            isRef.push_back(false);
    }
    // if (G->config->verbose)
    //     fprintf(stderr, "processing alleles for file %u\n", f);

    while (true) {
        nodeVal = G->get_node_last_char(nodeId);
        if (nodeVal == 0) {
            nodeId = G->fwd(nodeId);
            alleleSeqPos++;
            currStart++;
            continue;
        }
        seqVal = G->encode(seq.s[seqPos]);

        assert(nodeVal % G->alph_size == 6 || nodeVal % G->alph_size == seqVal);
        if (seqPos + 1 == seq.l)
            break;
        if (nodeVal % G->alph_size != 6) {
            seqPos++;
        } else {
            currStart++;
        }
        seqValNext = G->encode(seq.s[seqPos]) + 1;

        // find edge to next node
        out = G->outdegree(nodeId);
        //fprintf(stderr, "out %lu nodeId %lu\n", out, nodeId);
        if (out == 1) {
            if (branchMap.find(std::make_pair(nodeId, out)) != branchMap.end()) {
                currJoin = joins.at(branchMap[std::make_pair(nodeId, out)]);
                assert(currJoin.seqId1 == seqId);
                assert(currJoin.seqPos1 == alleleSeqPos);
                if (!isRefRun) {
                    if (alleleSeqPos - currStart + 1 > 0) {
                        SQLstream << "INSERT INTO AllelePathItem VALUES ("
                                  << f + 1 << ", " << alleleCnt++
                                           << ", " << seqId
                                           << ", " << currStart
                                           << ", " << alleleSeqPos - currStart + 1
                                           << ", 'TRUE')" << std::endl;
                    }
                } else {
                    isRef[seqId - 1] = true;                 
                }
                seqId = currJoin.seqId2;
                currStart = alleleSeqPos = currJoin.seqPos2;
                nodeId = G->outgoing(nodeId, seqValNext);
            } else {
                nodeId = G->fwd(nodeId);
                alleleSeqPos++;
            }
        } else {
            // find edge corresponding to the proceding seqPos
            uint64_t start = G->pred_last(nodeId - 1);
            uint64_t stop = G->succ_last(nodeId);
            assert(stop - start == out);
            edge = 0;
            size_t k;
            for (k = start + 1; k <= stop; ++k) {
                if (G->get_W(k) % G->alph_size == seqValNext) {
                    edge = k - start;
                    break;
                }
            }
            if (nodeVal % G->alph_size == 6)
                edge--;

            assert(edge > 0);
            if (branchMap.find(std::make_pair(nodeId, edge)) != branchMap.end()) {
                currJoin = joins.at(branchMap[std::make_pair(nodeId, edge)]);
                assert(currJoin.seqId1 == seqId);
                assert(currJoin.seqPos1 == alleleSeqPos);
                if (!isRefRun) {
                    if (alleleSeqPos - currStart + 1 > 0) {
                        SQLstream << "INSERT INTO AllelePathItem VALUES (" << f + 1
                                        << ", " << alleleCnt++
                                        << ", " << seqId
                                        << ", " << currStart
                                        << ", " << alleleSeqPos - currStart + 1
                                        << ", 'TRUE')" << std::endl;
                    }
                } else {
                    isRef[seqId - 1] = true;
                }
                seqId = currJoin.seqId2;
                currStart = alleleSeqPos = currJoin.seqPos2;
                nodeId = G->outgoing(nodeId, seqValNext);
            } else {
                nodeId = G->fwd(k);
                alleleSeqPos++;
            }
        }
    }
    if (!isRefRun) {
        if (alleleSeqPos - currStart + 1 > 0) {
            SQLstream << "INSERT INTO AllelePathItem VALUES (" << f + 1
                            << ", " << alleleCnt++
                            << ", " << seqId
                            << ", " << currStart
                            << ", " << alleleSeqPos - currStart + 1
                            << ", 'TRUE')" << std::endl;
        }
    } else {
        isRef[seqId - 1] = true;
        // write reference information to SQL stream 
        for (size_t i = 0; i < isRef.size(); ++i) {
            if (isRef[i])
                SQLstream << "INSERT INTO Reference VALUES ("
                  << i + 1 << ", 'seq"
                  << i + 1 << "', date('now'), "
                  << i + 1 << ", 0, NULL, NULL, NULL, NULL, NULL, 'TRUE');"
                  << std::endl;
            else
                SQLstream << "INSERT INTO Reference VALUES ("
                  << i + 1 << ", 'seq"
                  << i + 1 << "', date('now'), "
                  << i + 1 << ", 0, NULL, NULL, NULL, NULL, NULL, 'FALSE');"
                  << std::endl;
        }
        SQLstream << "INSERT INTO ReferenceSet VALUES (1, NULL, NULL, 'normal', 'FALSE');"
                  << std::endl;
        for (size_t i = 0; i < isRef.size(); ++i)
            SQLstream << "INSERT INTO Reference_ReferenceSet_Join VALUES ("
                            << i + 1 << ", 1);" << std::endl;
        for (size_t i = 0; i < joins.size(); ++i)
            SQLstream << "INSERT INTO GraphJoin_ReferenceSet_Join VALUES ("
                            << i + 1 << ", 1);" << std::endl;
    }

}


void toSQL(DBG_succ *G, const std::vector<std::string> &fname,
                        const std::string &sqlfbase) {
    // this vector stores the joins between the sequence objects we wrote
    std::vector<JoinInfo> joins;
    // we also store for each branching edge the join it creates
    // we will use this for allele traversal
    std::map<std::pair<uint64_t, DBG_succ::TAlphabet>, uint64_t> branchMap;

    // open sql filestream
    std::ofstream SQLstream(sqlfbase + ".sql");

    // traverse the graph, thereby filling joins vector, branchMap and 
    // writing the sequences to individual fasta files
    size_t seqNum = traverseGraph(G, joins, branchMap, sqlfbase, SQLstream); 
    
    // write graph joins to SQL file
    for (size_t i = 0; i < joins.size(); ++i) {
        if (joins.at(i).seqId1 < joins.at(i).seqId2 ||
             (joins.at(i).seqId1 == joins.at(i).seqId2 &&
              joins.at(i).seqPos1 < joins.at(i).seqPos2)) {
            SQLstream << "INSERT INTO GraphJoin VALUES ("
                      << i + 1 << ", " << joins.at(i).seqId1
                               << ", " << joins.at(i).seqPos1
                               << ", 'FALSE', " << joins.at(i).seqId2
                               << ", " << joins.at(i).seqPos2
                               << ", 'TRUE');" << std::endl;
        } else {
            SQLstream << "INSERT INTO GraphJoin VALUES ("
                      << i + 1 << ", " << joins.at(i).seqId2
                               << ", " << joins.at(i).seqPos2
                               << ", 'TRUE', " << joins.at(i).seqId1
                               << ", " << joins.at(i).seqPos1
                               << ", 'FALSE');" << std::endl;
        }
    }

    // for each input sequence traverse the graph once more and
    // collect allele path information 
    for (unsigned int f = 0; f < fname.size(); ++f) {

        // first traversal is for reference info
        if (f == 0) {
            // open stream to fasta file
            gzFile input_p = gzopen(fname.at(f).c_str(), "r");
            kseq_t *stream = kseq_init(input_p);

            if (stream)
                std::cerr << "ERROR while opening input file "
                          << fname.at(f) << std::endl;

            while (kseq_read(stream) >= 0) {
                allelesFromSeq(G, stream->seq, f, joins, branchMap, SQLstream, true, seqNum);
            }
            kseq_destroy(stream);
            gzclose(input_p);

            // open variant set
            SQLstream << "INSERT INTO VariantSet VALUES (1, 1, 'deBruijnGraph');"
                      << std::endl;
        }
        // open stream to fasta file
        gzFile input_p = gzopen(fname.at(f).c_str(), "r");
        kseq_t *stream = kseq_init(input_p);
        if (stream)
            std::cerr << "ERROR while opening input file "
                      << fname.at(f) << std::endl;

        while (kseq_read(stream) >= 0) 
            allelesFromSeq(G, stream->seq, f, joins, branchMap, SQLstream);

        kseq_destroy(stream);
        gzclose(input_p);
    }

    // write call set (one per input file)
    for (unsigned int f = 0; f < fname.size(); ++f)
        SQLstream << "INSERT INTO CallSet VALUES ("
                  << f + 1 << ", '" << fname.at(f)
                           << "', 'DBG"
                  << f + 1 << "');" << std::endl;
    for (unsigned int f = 0; f < fname.size(); ++f)
        SQLstream << "INSERT INTO VariantSet_CallSet_Join VALUES (1, "
                  << f + 1 << ");" << std::endl;
    for (unsigned int f = 0; f < fname.size(); ++f) {
        for (unsigned int ff = 0; ff < fname.size(); ++ff) {
            if (f == ff) {
                SQLstream << "INSERT INTO AlleleCall VALUES ("
                          << ff + 1 << ", " << f + 1 << ", 1);" << std::endl;
            } else {
                SQLstream << "INSERT INTO AlleleCall VALUES ("
                          << ff + 1 << ", " << f + 1 << ", 0);" << std::endl;
            }
        }
    }
}

}
