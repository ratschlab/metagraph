#include "dbg_succinct_merge.hpp"

#include <unistd.h>

#include "utils.hpp"
#include "dbg_succinct.hpp"


namespace merge {

/**
 * This object collects information about branches during graph traversal for the
 * purpose of merging, so we know where to jump back to when we reached a dead end.
 */
struct BranchInfoMerge {
    uint64_t nodeId;
    TAlphabet lastEdge;
    std::deque<TAlphabet> last_k;
};


BranchInfoMerge pop_branch(std::stack<BranchInfoMerge> &branchnodes,
                           uint64_t &nodeId, uint64_t &lastEdge,
                           std::deque<TAlphabet> &last_k);


void merge(DBG_succ *Gt,
           std::vector<DBG_succ*> Gv,
           std::vector<uint64_t> kv,
           std::vector<uint64_t> nv) {

    // Preliminarities
    for (size_t i = 0; i < Gv.size(); i++) {
        // check whether we can merge the given graphs
        if (i > 0 && (Gv.at(i)->get_k() != Gv.at(i-1)->get_k())) {
            fprintf(stderr, "Graphs have different k-mer lengths - cannot be merged!\n");
            exit(1);
        }
        // positions in the graph for respective traversal
        nv.at(i) = (nv.at(i) == 0) ? Gv.at(i)->get_W().size() : nv.at(i);
        // handle special cases where one or both input graphs are empty
        kv.at(i) = (kv.at(i) == 0) ? Gv.at(i)->get_W().size() : kv.at(i);
        //std::cerr << "k(" << i << ") " << kv.at(i) << " n(" << i << ") " << nv.at(i) << std::endl;
    }

    // keep track of how many nodes we added
    uint64_t added = 0;
    size_t cnt = 0;
    std::map<uint64_t, std::deque<TAlphabet>> last_added_nodes;
    // init last added nodes, if not starting from the beginning
    std::deque<TAlphabet> curr_seq;
    for (size_t i = 0; i < Gv.size(); i++) {
        if (kv.at(i) < 2)
            continue;
        for (size_t a = 0; a < Gt->alph_size; a++) {
            uint64_t sl = std::max(
                Gv.at(i)->pred_W(kv.at(i) - 1, a),
                Gv.at(i)->pred_W(kv.at(i) - 1, a + Gt->alph_size)
            );
            if (sl == 0)
                continue;
            auto la = last_added_nodes.find(a);
            curr_seq = Gv.at(i)->get_node_seq(sl);
            if (la == last_added_nodes.end()
                 || utils::colexicographically_greater(curr_seq, la->second))
                last_added_nodes[a] = curr_seq;
        }
    }

    if (Gt->verbose) {
        std::cout << "Size of bins to merge: " << std::endl;
        for (size_t i = 0; i < Gv.size(); i++) {
            std::cout << nv.at(i) - kv.at(i) << std::endl;
        }
    }

    // Send parallel pointers running through each of the graphs. At each step, compare all
    // graph nodes at the respective positions with each other. Insert the lexicographically
    // smallest one into the common merge graph G (this).
    while (true) {

        if (Gt->verbose && added > 0 && added % 1000 == 0) {
            std::cout << "." << std::flush;
            if (added % 10'000 == 0) {
                std::cout << "added " << added;
                for (size_t i = 0; i < Gv.size(); i++)
                    std::cout << " - G" << i << ": edge " << kv.at(i) << "/" << Gv.at(i)->get_W().size();
                std::cout << std::endl;
            }
        }

        // find set of smallest pointers
        std::pair<std::vector<bool>, uint64_t> smallest = utils::compare_nodes(Gv, kv, nv, cnt);
        if (cnt == 0)
            break;
        size_t curr_k = std::max_element(smallest.first.begin(), smallest.first.end())
                                 - smallest.first.begin();
        std::deque<TAlphabet> seq1 = Gv.at(curr_k)->get_node_seq(kv.at(curr_k));
        uint64_t val = Gv.at(curr_k)->get_W(kv.at(curr_k)) % Gt->alph_size;
        //std::cerr << "val: " << val << std::endl;
        //usleep(1000000);

        assert(val == smallest.second % Gt->alph_size);

        //std::cerr << "inserting into W" << std::endl;
        // check whether we already added a node whose outgoing edge points to the
        // same node as the current one
        auto it = last_added_nodes.find(smallest.second % Gt->alph_size);
        if (it != last_added_nodes.end() && utils::seq_equal(seq1, it->second, 1)) {
            Gt->W->insert(Gt->W->size(), val + Gt->alph_size);
        } else {
            Gt->W->insert(Gt->W->size(), smallest.second);
        }
        last_added_nodes[val] = seq1;
        Gt->update_F(Gv.at(curr_k)->get_node_last_char(kv.at(curr_k)), +1);
        Gt->last->insertBit(Gt->W->size() - 1, true);

        // handle multiple outgoing edges
        if (added > 0 && Gt->W->size() > 2 && (*(Gt->last))[Gt->W->size() - 2]) {
            // compare the last two added nodes
            auto it1 = last_added_nodes.find((*(Gt->W))[Gt->W->size() - 2] % Gt->alph_size);
            auto it2 = last_added_nodes.find((*(Gt->W))[Gt->W->size() - 1] % Gt->alph_size);
            if (it1 != last_added_nodes.end()
                     && it2 != last_added_nodes.end()
                     && it1 != it2
                     && utils::seq_equal(it1->second, it2->second)) {
                Gt->last->set(Gt->W->size() - 2, false);
            }
        }
        uint64_t updated = 0;
        for (size_t i = 0; i < Gv.size(); i++) {
            if (smallest.first.at(i)) {
                updated += (kv.at(i) < nv.at(i));
                kv.at(i) += (kv.at(i) < nv.at(i));
            }
        }
        ++added;
        if (updated == 0)
            break;
    }
    Gt->p_ = Gt->succ_W(1, 0);
}

/**
* Heavily borrowing from the graph sequence traversal, this function gets a graph pointer Gm and merges its
* nodes into the target graph object Gt. The edges of Gm are fully traversed and nodes are added to
* Gt if not existing yet. This function is well suited to merge small graphs into large ones.
*/
void merge(DBG_succ *Gt, DBG_succ *Gm) {

    // store all branch nodes on the way
    std::stack<BranchInfoMerge> branchnodes;
    // bool vector that keeps track of visited nodes
    std::vector<bool> visited(Gm->get_W().size(), false);

    // some initializations
    uint64_t nodeId = 1; // start at source node
    size_t out = Gm->outdegree(nodeId);
    BranchInfoMerge branch;
    TAlphabet val;
    TAlphabet lastEdge = 0;
    // keep a running list of the last k-1 characters we have seen
    std::deque<TAlphabet> last_k;
    bool new_branch = false;
    bool old_last = (*Gt->last)[Gt->p_];
    bool initial_k = true;
    uint64_t added = 0;

    // keep traversing until we reach the think and have worked off all branches from the stack
    while (out > 0 || branchnodes.size() > 0) {

        // verbose output
        if (added > 0 && added % 1000 == 0) {
            std::cout << "." << std::flush;
            if (added % 10000 == 0) {
                fprintf(stdout, "merged %" PRIu64
                                    " / %" PRIu64
                              " - edges %" PRIu64
                              " / nodes %" PRIu64 "\n",
                                added,
                                Gm->get_W().size(),
                                Gt->W->size() - 1,
                                Gt->rank_last((Gt->last->size() - 1)));
            }
        }

        // we have reached the sink but there are unvisited nodes left on the stack
        if (out == 0) {
            if (branchnodes.size() == 0)
                break;
            // get new branch
            branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
            out = Gm->outdegree(nodeId);
            new_branch = true;
        }
        //std::cerr << "starting loop with nodID " << nodeId << std::endl;

        if (new_branch) {
            // find node where to restart insertion
            uint64_t ridx = Gt->index(last_k, last_k.size());
            // put at the beginning of equal node range
            ridx = Gt->pred_last(ridx - 1) + 1;
            ridx -= (Gt->p_ < ridx);
            //std::cerr << "ridx: " << ridx << std::endl;
            //std::cerr << "bef move " << std::endl;
            //this->print_seq();
            assert(!(*(Gt->last))[Gt->p_]); // must not remove a dangling end
            // move p to last position
            Gt->W->remove(Gt->p_);
            Gt->last->deleteBit(Gt->p_);
            Gt->update_F(Gt->get_node_last_char(Gt->p_), -1);
            //std::cerr << "moving p from " << p << " to " << ridx << std::endl;
            Gt->p_ = ridx;
            Gt->W->insert(Gt->p_, 0);
            Gt->last->insertBit(Gt->p_, 0); // it's at the beginning of a branch node
            Gt->update_F(Gt->get_node_last_char(Gt->p_), +1);

            new_branch = false;
            //std::cerr << "aft move " << std::endl;
            //this->print_seq();
        }

        // we have not visited that node before
        if (!visited.at(nodeId)) {
            visited.at(nodeId) = true;
            val = Gm->get_W(nodeId) % Gt->alph_size;
            //std::cerr << "current val " << val % alph_size
            //          << " nodeID: " << nodeId << std::endl;
            //G->print_seq();
            last_k.push_back(Gm->get_node_last_char(nodeId));
            if (last_k.size() > Gt->k_)
                last_k.pop_front();
        }

        // there is only one child
        if (out == 1) {
            uint64_t next = Gm->fwd(nodeId);
            val = Gm->get_W(nodeId) % Gt->alph_size;
            if ((val != 6 || !initial_k) && val != 0) {
                initial_k = false;
                Gt->append_pos(val % Gt->alph_size);
                added++;
                //std::cerr << "append " << val % alph_size << " nodeID: " << nodeId << std::endl;
                //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                //this->print_seq();
            }

            // the next node is new
            if (!visited.at(next)) {
                nodeId = next;
                lastEdge = 0;
            // we have seen the next node before
            } else {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // append next node
                if ((*(Gt->last))[Gt->p_]) {
                    val = Gm->get_W(next) % Gt->alph_size;
                    if ((val != 6 || !initial_k) && val != 0) {
                        initial_k = false;
                        Gt->append_pos(val % Gt->alph_size);
                        added++;
                        //std::cerr << "..append " << val % alph_size << " nodeID: " << nodeId << std::endl;
                        //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                    }
                }
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = Gm->outdegree(nodeId);
                //std::cerr << "new branch 1 - nodeID: " << nodeId << " next is " << next << std::endl;
                new_branch = true;
            }
        // there are several children
        } else {
            size_t cnt = 0;
            bool updated = false;
            // account for sentinel symbol as possible outgoing edge
            if (Gm->get_W(nodeId) == 0)
                cnt++;
            // loop over outgoing edges
            for (TAlphabet c = 1; c < Gt->alph_size; ++c) {
                uint64_t next = Gm->outgoing(nodeId, c);
                if (next > 0) {
                    cnt++;
                    // we already handled this edge erlier
                    if (cnt <= lastEdge)
                        continue;
                    lastEdge++;

                    //std::cerr << "mult edge - val is now " << c << std::endl;
                    uint64_t curr_p = Gt->p_;
                    if ((c != 6 || !initial_k) && c != 0) {
                        initial_k = false;
                        //std::cerr << "p (before): " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                        Gt->append_pos(c);
                        added++;
                        curr_p += (Gt->p_ <= curr_p);
                        //std::cerr << "append " << c % alph_size << " nodeID: " << nodeId << std::endl;
                        //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                    }

                    if (!visited.at(next)) {
                        // there are remaining branches - push node to stack
                        if (cnt < out && next != nodeId) {
                            // push node information to stack
                            branchnodes.push({ nodeId, lastEdge, last_k });
                            //std::cerr << "pushing nodeID " << nodeId << " onto stack" << std::endl;
                        }
                        nodeId = next;
                        updated = true;
                        lastEdge = 0;
                        break;
                    } else {
                        //std::cerr << "visited next before: " << next <<std::endl;
                        // append next node
                        if ((*(Gt->last))[Gt->p_]) {
                            c = Gm->get_W(next) % Gt->alph_size;
                            if ((c != 6 || !initial_k) && c != 0) {
                                initial_k = false;
                                Gt->append_pos(c);
                                added++;
                                //std::cerr << "...append " << c % alph_size << " nodeID: " << nodeId << std::endl;
                                //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                                //this->print_seq();
                            }
                        }
                        // reset to previous position
                        if (nodeId != next) {
                            Gt->W->remove(Gt->p_);
                            Gt->last->deleteBit(Gt->p_);
                            Gt->update_F(Gt->get_node_last_char(Gt->p_), -1);
                            curr_p -= (Gt->p_ < curr_p);
                            //std::cerr << ".moving p from " << p << " to " << curr_p << std::endl;
                            Gt->p_ = curr_p;
                            Gt->W->insert(Gt->p_, 0);
                            Gt->last->insertBit(Gt->p_, 0); // it's at the beginning of a branch node
                            Gt->update_F(Gt->get_node_last_char(Gt->p_), +1);
                            //std::cerr << ".aft move " << std::endl;
                            //this->print_seq();
                        }
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
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = Gm->outdegree(nodeId);
                new_branch = true;
                //std::cerr << "new branch 2" << std::endl;
            }
        }
        out = Gm->outdegree(nodeId);
    }
    // bring graph into default state
    std::deque<TAlphabet> tmp_p;
    for (size_t t = 0; t < Gt->k_; t++) {
        tmp_p.push_back(6);
    }
    uint64_t old_p = Gt->pred_last(Gt->index(tmp_p, tmp_p.size()) - 1) + 1;

    old_p -= (Gt->p_ < old_p);
    Gt->W->remove(Gt->p_);
    Gt->last->deleteBit(Gt->p_);
    Gt->update_F(Gt->get_node_last_char(Gt->p_), -1);
    Gt->p_ = old_p;
    Gt->W->insert(Gt->p_, 0);
    Gt->last->insertBit(Gt->p_, old_last);

    // locally update sorting
    std::pair<uint64_t, uint64_t> R = Gt->get_equal_node_range(Gt->p_);
    if (R.second - R.first > 0) {
        Gt->sort_W_locally(R.first, R.second);
        while ((*(Gt->W))[Gt->p_] != 0) {
            (Gt->p_)--;
        }
        assert((*(Gt->W))[Gt->p_] == 0);
    }
    Gt->update_F(Gt->get_node_last_char(Gt->p_), +1);
}


BranchInfoMerge pop_branch(std::stack<BranchInfoMerge> &branchnodes,
                           uint64_t &nodeId, uint64_t &lastEdge,
                           std::deque<TAlphabet> &last_k) {
    BranchInfoMerge branch = branchnodes.top();
    branchnodes.pop();
    lastEdge = branch.lastEdge;
    nodeId = branch.nodeId;
    last_k = branch.last_k;

    return branch;
}


void traversalHash(DBG_succ *G) {

    // store all branch nodes on the way
    std::stack<BranchInfoMerge> branchnodes;
    // bool vector that keeps track of visited nodes
    std::vector<bool> visited(G->get_W().size(), false);

    // some initializations
    uint64_t nodeId = 1; // start at source node
    uint64_t count = 0;
    size_t out = G->outdegree(nodeId);
    BranchInfoMerge branch;
    TAlphabet lastEdge = 0;
    // keep a running list of the last k-1 characters we have seen
    std::deque<TAlphabet> last_k;

    // keep traversing until we reach the sink and have worked off all branches from the stack
    while (out > 0 || branchnodes.size() > 0) {

        if (count > 0 && count % 100'000 == 0) {
            std::cout << "." << std::flush;
            if (count % 1'000'000 == 0)
                std::cout << count << std::endl;
        }

        // we have reached the sink but there are unvisited nodes left on the stack
        if (out == 0) {
            if (branchnodes.size() == 0)
                break;

            // get new branch
            branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
            out = G->outdegree(nodeId);
        }

        // we have not visited that node before
        if (!visited.at(nodeId)) {
            visited.at(nodeId) = true;
            last_k.push_back(G->get_node_last_char(nodeId));
            if (last_k.size() < G->k_)
                last_k = G->get_node_seq(nodeId);
            if (last_k.size() > G->k_)
                last_k.pop_front();
        }

        // there is only one child
        if (out == 1) {

            //std::cout << " " << get_alphabet_symbol(this->get_W(nodeId) % alph_size) << " " << nodeId << std::endl;
            count++;
            uint64_t next = G->fwd(nodeId);
            // the next node is new
            if (!visited.at(next)) {
                nodeId = next;
                lastEdge = 0;
            // we have seen the next node before
            // --> jump back to previous branch
            } else {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = G->outdegree(nodeId);
            }
        // there are several children
        } else {
            // account for sentinel symbol as possible outgoing edge
            size_t cnt = (G->get_W(nodeId) == 0);
            bool updated = false;
            // loop over outgoing edges
            for (TAlphabet c = 1; c < G->alph_size; ++c) {
                uint64_t next = G->outgoing(nodeId, c);
                if (next > 0) {
                    cnt++;
                    // we already handled this edge erlier
                    if (cnt <= lastEdge)
                        continue;
                    lastEdge++;

                    //std::cout << " " << get_alphabet_symbol(c) << " " << pred_W(nodeId, c) << std::endl;
                    count++;

                    if (!visited.at(next)) {
                        // there are remaining branches - push node to stack
                        if (cnt < out && next != nodeId) {
                            // push node information to stack
                            branchnodes.push({ nodeId, lastEdge, last_k });
                        }
                        nodeId = next;
                        updated = true;
                        lastEdge = 0;
                        break;
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
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = G->outdegree(nodeId);
            }
        }
        out = G->outdegree(nodeId);
    }
    // handle current end
    std::cout << " " << G->decode(0) << " " << G->p_;
    std::cout << std::endl;
}

} // namespace merge
