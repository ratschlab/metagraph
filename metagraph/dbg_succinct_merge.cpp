#include "dbg_succinct_merge.hpp"

#include <unistd.h>
#include <mutex>

#include "utils.hpp"
#include "dbg_succinct.hpp"
#include "config.hpp"


namespace merge {


/**
 * Helper function to determine the bin boundaries, given
 * a number of bins basen on the total number of nodes in graph.
 * [v[0], v[1]), [v[1], v[2]), ...
 */
std::vector<uint64_t> get_bins(const DBG_succ *G, uint64_t num_bins) {
    assert(num_bins > 0);

    uint64_t num_nodes = G->rank_last(G->get_W().size() - 1);

    if (G->verbose && num_bins > num_nodes) {
        std::cerr << "WARNING: There are max "
                  << num_nodes << " slots available for binning. Your current choice is "
                  << num_bins << " which will create "
                  << num_bins - num_nodes << " empty slots." << std::endl;
    }

    std::vector<uint64_t> result { 1 };

    uint64_t binsize = num_nodes / num_bins + (num_nodes % num_bins != 0);
    uint64_t num_full_bins = num_nodes / binsize;

    for (uint64_t i = 1; i <= num_full_bins; ++i) {
        result.push_back(G->select_last(i * binsize) + 1);
    }
    if (num_nodes % num_bins) {
        result.push_back(G->get_W().size());
    }

    while (result.size() < num_bins + 1) {
        result.push_back(result.back());
    }

    return result;
}

/**
 * Subset the bins to the chunk
 * computed in the current distributed compute.
 */
std::vector<uint64_t> subset_bins(const std::vector<uint64_t> &ref_bins,
                                  size_t chunk_idx, size_t num_chunks) {
    assert(ref_bins.size() % num_chunks == 0);

    size_t begin = chunk_idx * (ref_bins.size() / num_chunks);
    size_t end = (chunk_idx + 1) * (ref_bins.size() / num_chunks);

    return std::vector<uint64_t>(ref_bins.begin() + begin,
                                 ref_bins.begin() + end);
}

std::vector<uint64_t> get_chunk(const DBG_succ &graph,
                                const std::vector<std::deque<TAlphabet>> &border_kmers,
                                bool with_tail) {
    std::vector<uint64_t> result;

    for (size_t i = 0; i < border_kmers.size(); i++) {

        uint64_t last = graph.pred_kmer(border_kmers[i]);

        if (graph.get_node_seq(last) != border_kmers[i]) {
            result.push_back(last + 1);
        } else {
            result.push_back(graph.pred_last(last - 1) + 1);
        }
    }
    if (with_tail) {
        result.push_back(graph.get_W().size());
    }
    return result;
}

/**
 * Show an overview of the distribution of merging bin sizes.
 */
void print_bin_stats(const std::vector<std::vector<uint64_t>> &bins) {
    size_t min_bin = 0;
    size_t max_bin = 0;
    size_t total_bin = 0;
    size_t num_bins = 0;

    for (size_t j = 0; j < bins.size(); j++) {
        num_bins = bins[j].size() - 1;
        size_t cum_size = 0;
        for (size_t i = 0; i < bins.front().size() - 1; ++i) {
            cum_size += bins.at(j).at(i + 1) - bins.at(j).at(i);
        }
        if (cum_size > 0) {
            min_bin = (min_bin == 0)
                        ? cum_size
                        : std::min(min_bin, cum_size);
            max_bin = std::max(max_bin, cum_size);
        }
        total_bin += cum_size;
    }

    std::cout << "\nTotal number of bins: " << num_bins << "\n";
    if (num_bins) {
        std::cout << "Total size: " << total_bin << "\n";
        std::cout << "Smallest bin: " << min_bin << "\n";
        std::cout << "Largest bin: " << max_bin << "\n";
        std::cout << "Average bin size: " << total_bin / num_bins << "\n";
    }
    std::cout << std::endl;
}

struct ParallelMergeContext {
    size_t idx = 0;

    std::mutex mu;
};

/**
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used.
 */
void parallel_merge_wrapper(const std::vector<const DBG_succ*> &graphs,
                            const std::vector<std::vector<uint64_t>> &bins,
                            ParallelMergeContext *context,
                            std::vector<const DBG_succ*> *merged_chunks) {
    assert(merged_chunks && context);
    assert(graphs.size() > 0);
    assert(graphs.size() == bins.size());
    assert(merged_chunks->size() == bins.front().size() - 1);

    while (true) {
        context->mu.lock();

        if (context->idx == bins.front().size() - 1) {
            context->mu.unlock();
            break;
        }

        size_t curr_idx = context->idx++;

        context->mu.unlock();

        std::vector<uint64_t> kv;
        std::vector<uint64_t> nv;
        for (size_t i = 0; i < graphs.size(); i++) {
            kv.push_back(bins.at(i).at(curr_idx));
            nv.push_back(bins.at(i).at(curr_idx + 1));
        }
        DBG_succ *merged_chunk = merge::merge(graphs, kv, nv);
        assert(merged_chunk && "Merge gone wrong");

        merged_chunks->at(curr_idx) = merged_chunk;
    }
}

// TODO: return plain vectors instead of succinct graph
//       because it has to be merged with other chunks afterwards anyway.
DBG_succ* build_chunk(const std::vector<const DBG_succ*> &graphs,
                      size_t chunk_idx,
                      size_t num_chunks,
                      size_t num_threads,
                      size_t num_bins_per_thread) {
    assert(graphs.size() > 0);
    assert(num_chunks > 0);
    assert(chunk_idx < num_chunks);

    // get bins in graphs according to required threads
    if (graphs.front()->verbose) {
        std::cout << "Collecting reference bins" << std::endl;
        std::cout << "parallel " << num_threads
                  << " per thread " << num_bins_per_thread
                  << " parts total " << num_chunks << std::endl;
    }

    auto ref_bins = subset_bins(
        get_bins(graphs.front(),
                 num_threads * num_bins_per_thread * num_chunks),
        chunk_idx,
        num_chunks
    );
    std::vector<std::deque<TAlphabet>> border_kmers;
    bool with_tail = false;
    for (size_t i = 0; i < ref_bins.size(); ++i) {
        if (ref_bins[i] == graphs.front()->get_W().size()) {
            with_tail = true;
            break;
        }
        border_kmers.push_back(graphs.front()->get_node_seq(ref_bins[i]));
        // Make the kmers from different bins differ
        // by a prefix of length at least 2.
        // So that we don't have to adjust the W array when concatenating.
        border_kmers.back()[0] = DBG_succ::encode('$');
    }

    if (graphs.front()->verbose)
        std::cout << "Collecting relative bins" << std::endl;

    std::vector<std::vector<uint64_t>> bins;
    for (size_t i = 0; i < graphs.size(); i++) {
        bins.push_back(get_chunk(*graphs.at(i), border_kmers, with_tail));
    }

    // print bin stats
    if (graphs.front()->verbose)
        print_bin_stats(bins);

    // create threads and start the jobs
    std::vector<std::thread> threads;
    ParallelMergeContext context;
    std::vector<const DBG_succ*> merged_chunks(bins.front().size() - 1, NULL);
    for (size_t tid = 0; tid < num_threads; tid++) {
        threads.emplace_back(parallel_merge_wrapper, graphs,
                                                     bins,
                                                     &context,
                                                     &merged_chunks);
        if (graphs.front()->verbose)
            std::cout << "starting thread " << tid << std::endl;
    }

    // join threads
    if (graphs.front()->verbose)
        std::cout << "Waiting for threads to join" << std::endl;

    for (size_t tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }

    // collect results
    if (graphs.front()->verbose)
        std::cout << "Collecting results" << std::endl;

    DBG_succ *graph = new DBG_succ(graphs.front()->get_k());
    // remove dummy source
    graph->W->remove(1);
    graph->last->deleteBit(1);
    graph->update_F(DBG_succ::encode('$'), -1);
    uint64_t cur_pos = 1;

    for (const DBG_succ *merged_chunk : merged_chunks) {
        assert(merged_chunk);

        graph->verbose_cout("    adding ", merged_chunk->W->size(), " edges\n");

        // handle last and W
        for (uint64_t j = 1; j < merged_chunk->W->size(); ++j) {
            graph->W->insert(cur_pos, merged_chunk->get_W(j));
            graph->last->insertBit(cur_pos, merged_chunk->get_last(j));
            cur_pos++;
        }

        graph->verbose_cout("new total edges: ", graph->W->size(), "\n");

        // handle F
        assert(graph->F.size() == merged_chunk->F.size());

        for (size_t p = 0; p < graph->F.size(); ++p) {
            graph->F.at(p) += merged_chunk->F[p];
        }

        delete merged_chunk;
    }

    return graph;
}

// TODO: save/load plain vectors instead of succinct graphs
DBG_succ* merge_chunks(const std::string &filenamebase, size_t num_chunks) {
    // the bit array indicating the last outgoing edge of a node (static container for full init)
    std::vector<bool> last_stat { 0 };

    // the array containing the edge labels
    std::vector<TAlphabet> W_stat { 0 };

    std::vector<uint64_t> F(DBG_succ::alph_size, 0);

    size_t k = 0;

    for (uint64_t f = 0; f < num_chunks; f++) {
        std::string filename = filenamebase
                                + "." + std::to_string(f)
                                + "_" + std::to_string(num_chunks);
        std::cout << "Opening file " << filename << std::endl;

        DBG_succ graph_to_append;
        if (!graph_to_append.load(filename)) {
            std::cerr << "ERROR: input file "
                      << filename << " corrupted" << std::endl;
            return NULL;
        }
        k = graph_to_append.get_k();

        graph_to_append.verbose_cout("    adding ", graph_to_append.W->size(), " edges\n");

        assert(dynamic_cast<wavelet_tree_dyn*>(graph_to_append.W));
        auto G_W_stat = dynamic_cast<wavelet_tree_dyn*>(graph_to_append.W)->to_vector();

        W_stat.insert(W_stat.end(), G_W_stat.begin() + 1, G_W_stat.end());

        for (size_t i = 1; i < graph_to_append.W->size(); ++i) {
            last_stat.push_back(graph_to_append.get_last(i));
        }

        graph_to_append.verbose_cout("new total edges: ", W_stat.size(), "\n");

        // handle F
        assert(F.size() == graph_to_append.F.size());
        for (size_t j = 0; j < F.size(); ++j) {
            F.at(j) += graph_to_append.F.at(j);
        }
    }

    DBG_succ *graph = new DBG_succ(k);

    delete graph->W;
    graph->W = new wavelet_tree_dyn(4, W_stat);

    delete graph->last;
    graph->last = new bit_vector_dyn(last_stat);

    graph->F = std::move(F);

    return graph;
}


DBG_succ* merge(const std::vector<const DBG_succ*> &Gv) {
    std::vector<uint64_t> kv;
    std::vector<uint64_t> nv;

    for (size_t i = 0; i < Gv.size(); ++i) {
        kv.push_back(1);
        nv.push_back(Gv[i]->get_W().size());
    }

    DBG_succ *merged = merge(Gv, kv, nv);
    return merged;
}


// TODO: return result in plain vectors, construct a DBG_succ from them if needed
DBG_succ* merge(const std::vector<const DBG_succ*> &Gv,
                std::vector<uint64_t> kv,
                std::vector<uint64_t> nv) {

    assert(kv.size() == Gv.size());
    assert(nv.size() == Gv.size());

    DBG_succ *Gt = new DBG_succ(Gv.at(0)->get_k());
    Gt->W->remove(1);
    Gt->last->deleteBit(1);
    Gt->update_F(DBG_succ::encode('$'), -1);

    for (size_t i = 0; i < Gv.size(); i++) {
        // check whether we can merge the given graphs
        assert(Gv.at(i)->get_k() == Gt->get_k()
                && "Graphs have different k-mer lengths - cannot be merged!\n");

        // positions in the graph for respective traversal
        nv.at(i) = (nv.at(i) == 0) ? Gv.at(i)->get_W().size() : nv.at(i);
        // handle special cases where one or both input graphs are empty
        kv.at(i) = (kv.at(i) == 0) ? Gv.at(i)->get_W().size() : kv.at(i);
        //std::cerr << "k(" << i << ") " << kv.at(i)
        //          << " n(" << i << ") " << nv.at(i) << std::endl;
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
                    std::cout << " - G" << i << ": edge " << kv.at(i)
                                             << "/" << Gv.at(i)->get_W().size();
                std::cout << std::endl;
            }
        }

        // find set of smallest pointers
        auto smallest = utils::compare_nodes(Gv, kv, nv, cnt);

        if (cnt == 0)
            break;
        size_t curr_k = std::max_element(smallest.first.begin(), smallest.first.end())
                         - smallest.first.begin();

        auto seq1 = Gv.at(curr_k)->get_node_seq(kv.at(curr_k));

        uint64_t val = Gv.at(curr_k)->get_W(kv.at(curr_k)) % Gt->alph_size;
        //std::cerr << "val: " << val << std::endl;
        //usleep(1000000);

        assert(val == smallest.second % Gt->alph_size);

        //std::cerr << "inserting into W" << std::endl;
        // check whether we already added a node whose outgoing edge points to the
        // same node as the current one
        auto it = last_added_nodes.find(smallest.second % Gt->alph_size);
        if (it != last_added_nodes.end()
                        && utils::seq_equal(seq1, it->second, 1)
                        && val != DBG_succ::encode('$')) {
            Gt->W->insert(Gt->W->size(), val + Gt->alph_size);
        } else {
            Gt->W->insert(Gt->W->size(), smallest.second);
        }
        last_added_nodes[val] = seq1;
        Gt->update_F(Gv.at(curr_k)->get_node_last_char(kv.at(curr_k)), +1);
        Gt->last->insertBit(Gt->W->size() - 1, true);

        // handle multiple outgoing edges
        if (added > 0 && Gt->W->size() > 2 && Gt->get_last(Gt->W->size() - 2)) {
            // compare the last two added nodes
            auto it1 = last_added_nodes.find(Gt->get_W(Gt->W->size() - 2) % Gt->alph_size);
            auto it2 = last_added_nodes.find(Gt->get_W(Gt->W->size() - 1) % Gt->alph_size);
            if (it1 != last_added_nodes.end()
                     && it2 != last_added_nodes.end()
                     && it1 != it2
                     && utils::seq_equal(it1->second, it2->second)) {
                if (it2->second.back() != DBG_succ::encode('$')
                        && Gt->get_W(Gt->W->size() - 2) == DBG_succ::encode('$')) {
                    Gt->update_F(Gv.at(curr_k)->get_node_last_char(kv.at(curr_k)), -1);
                    Gt->last->deleteBit(Gt->W->size() - 2);
                    Gt->W->remove(Gt->W->size() - 2);
                } else {
                    Gt->last->set(Gt->W->size() - 2, false);
                }
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

    return Gt;
}


} // namespace merge
