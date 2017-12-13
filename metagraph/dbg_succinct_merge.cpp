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
    uint64_t node_id;
    std::deque<TAlphabet> source_kmer;
};


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
}


void pop_branch(std::stack<BranchInfoMerge> *branchnodes,
                uint64_t *node_id,
                std::deque<TAlphabet> *source_kmer) {
    BranchInfoMerge &branch = branchnodes->top();
    *node_id = branch.node_id;
    *source_kmer = branch.source_kmer;
    branchnodes->pop();
}


/**
* Heavily borrowing from the graph sequence traversal, this function gets a graph pointer Gm and merges its
* nodes into the target graph object Gt. The edges of Gm are fully traversed and nodes are added to
* Gt if not existing yet. This function is well suited to merge small graphs into large ones.
*/
void merge(DBG_succ *Gt, DBG_succ *Gm) {
    // bool vector that keeps track of visited nodes
    std::vector<bool> marked(Gm->rank_last(Gm->get_W().size() - 1), false);

    // start at the source node
    uint64_t Gt_source_node = 1;
    uint64_t Gm_source_node = 1;
    // keep a running list of the last k characters we have seen
    std::deque<TAlphabet> source_kmer(Gt->get_k(), DBG_succ::encode('$'));

    // store all branch nodes on the way
    std::stack<BranchInfoMerge> branchnodes;
    branchnodes.push({ Gm_source_node, source_kmer });
    marked[Gm_source_node] = true;
    bool new_branch_started = true;

    uint64_t added_counter = 0;

    // keep traversing until we reach the sink and have worked off all branches from the stack
    while (branchnodes.size()) {
        // get new branch
        pop_branch(&branchnodes, &Gm_source_node, &source_kmer);

        if (new_branch_started) {
            // find node where to restart insertion
            Gt_source_node = Gt->index(source_kmer, source_kmer.size());
            // put at the beginning of equal node range
            // Gt_source_node = Gt->pred_last(Gt_source_node - 1) + 1;
            new_branch_started = false;
        }
        // verbose output
        if (added_counter > 0 && added_counter % 1000 == 0) {
            std::cout << "." << std::flush;
            if (added_counter % 10000 == 0) {
                std::cout << "merged " << std::to_string(added_counter)
                          << " / " << std::to_string(Gm->get_W().size())
                          << " - edges " << std::to_string(Gt->W->size() - 1)
                          << " / nodes " << std::to_string(Gt->rank_last((Gt->last->size() - 1)))
                          << "\n";
            }
        }

        // we have reached the sink but there are unvisited nodes left on the stack
        if (Gm_source_node > 1 && Gm->get_W(Gm_source_node) == DBG_succ::encode('$')) {
            new_branch_started = true;
            continue;
        }

        std::deque<TAlphabet> target_kmer = source_kmer;
        target_kmer.pop_front();
        target_kmer.push_back(0);

        // loop over outgoing edges
        for (TAlphabet c = 1; c < Gt->alph_size; ++c) {
            uint64_t target_node = Gm->outgoing(Gm_source_node, c);
            if (!target_node)
                continue;

            //this->print_seq();
            Gt_source_node = Gt->append_pos(c, Gt_source_node);
            added_counter++;
            //std::cerr << "append " << c % alph_size
            //          << " Gm_source_node: " << Gm_source_node << std::endl;
            //this->print_seq();

            if (marked.at(target_node))
                continue;
            marked.at(target_node) = true;

            // push node information to stack
            target_kmer[Gt->get_k() - 1] = c;
            branchnodes.push({ target_node, target_kmer });

            new_branch_started = false;
        }
    }
}

} // namespace merge
