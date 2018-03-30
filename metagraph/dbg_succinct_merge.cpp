#include "dbg_succinct_merge.hpp"

#include <thread>
#include <mutex>

#include "utils.hpp"


namespace merge {

/**
 * Helper function to determine the bin boundaries, given
 * a number of bins basen on the total number of nodes in graph.
 * [v[0], v[1]), [v[1], v[2]), ...
 */
std::vector<uint64_t> get_bins(const DBG_succ &G, uint64_t num_bins, bool verbose) {
    assert(num_bins > 0);

    uint64_t num_nodes = G.num_nodes();

    if (verbose && num_bins > num_nodes) {
        std::cerr << "WARNING: There are max "
                  << num_nodes << " slots available for binning. Your current choice is "
                  << num_bins << " which will create "
                  << num_bins - num_nodes << " empty slots." << std::endl;
    }

    std::vector<uint64_t> result { 1 };

    uint64_t binsize = num_nodes / num_bins + (num_nodes % num_bins != 0);
    uint64_t num_full_bins = num_nodes / binsize;

    for (uint64_t i = 1; i <= num_full_bins; ++i) {
        result.push_back(G.select_last(i * binsize) + 1);
    }
    if (num_nodes % num_bins) {
        result.push_back(G.get_W().size());
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
    assert(ref_bins.size() > 1);
    assert((ref_bins.size() - 1) % num_chunks == 0);
    assert(chunk_idx < num_chunks);

    size_t chunk_size = (ref_bins.size() - 1) / num_chunks;

    return std::vector<uint64_t>(ref_bins.begin() + chunk_idx * chunk_size,
                                 ref_bins.begin() + (chunk_idx + 1) * chunk_size + 1);
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


void merge_blocks(const std::vector<const DBG_succ*> &Gv,
                  std::vector<uint64_t> kv,
                  const std::vector<uint64_t> &nv,
                  DBG_succ::Chunk *chunk,
                  bool verbose);

/**
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used.
 */
void parallel_merge_wrapper(const std::vector<const DBG_succ*> &graphs,
                            const std::vector<std::vector<uint64_t>> &bins,
                            ParallelMergeContext *context,
                            std::vector<DBG_succ::Chunk*> *chunks,
                            bool verbose) {
    assert(context);
    assert(graphs.size() > 0);
    assert(graphs.size() == bins.size());
    assert(chunks->size() == bins.front().size() - 1);

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
        merge_blocks(graphs, kv, nv, chunks->at(curr_idx), verbose);
    }
}


DBG_succ::Chunk* stack_graph_chunks(const std::vector<DBG_succ::Chunk*> &graph_chunks) {
    DBG_succ::Chunk *merged = new DBG_succ::Chunk();

    for (DBG_succ::Chunk *chunk : graph_chunks) {
        assert(chunk);
        merged->extend(*chunk);
        delete chunk;
    }

    return merged;
}


DBG_succ::Chunk* merge_blocks_to_chunk(const std::vector<const DBG_succ*> &graphs,
                                       size_t chunk_idx,
                                       size_t num_chunks,
                                       size_t num_threads,
                                       size_t num_bins_per_thread,
                                       bool verbose) {
    assert(graphs.size() > 0);
    assert(num_chunks > 0);
    assert(chunk_idx < num_chunks);

    // get bins in graphs according to required threads
    if (verbose) {
        std::cout << "Collecting reference bins" << std::endl;
        std::cout << "parallel " << num_threads
                  << " per thread " << num_bins_per_thread
                  << " parts total " << num_chunks << std::endl;
    }

    auto ref_bins = subset_bins(
        get_bins(*graphs.front(),
                 num_threads * num_bins_per_thread * num_chunks,
                 verbose),
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

    if (verbose)
        std::cout << "Collecting relative bins" << std::endl;

    std::vector<std::vector<uint64_t>> bins;
    for (size_t i = 0; i < graphs.size(); i++) {
        bins.push_back(get_chunk(*graphs.at(i), border_kmers, with_tail));
    }

    // print bin stats
    if (verbose)
        print_bin_stats(bins);

    // create threads and start the jobs
    std::vector<std::thread> threads;
    ParallelMergeContext context;

    std::vector<DBG_succ::Chunk*> blocks(bins.front().size() - 1);
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = new DBG_succ::Chunk();
    }

    for (size_t tid = 0; tid < num_threads; tid++) {
        threads.emplace_back(parallel_merge_wrapper, graphs,
                                                     bins,
                                                     &context,
                                                     &blocks,
                                                     verbose);
        if (verbose)
            std::cout << "starting thread " << tid << std::endl;
    }

    // join threads
    if (verbose)
        std::cout << "Waiting for threads to join" << std::endl;

    for (size_t tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }

    // collect results
    if (verbose)
        std::cout << "Collecting results" << std::endl;

    return stack_graph_chunks(blocks);
}


DBG_succ* merge(const std::vector<const DBG_succ*> &Gv, bool verbose) {
    std::vector<uint64_t> kv;
    std::vector<uint64_t> nv;

    for (size_t i = 0; i < Gv.size(); ++i) {
        kv.push_back(1);
        nv.push_back(Gv[i]->get_W().size());
    }

    auto *merged = new DBG_succ::Chunk();

    merge_blocks(Gv, kv, nv, merged, verbose);

    return DBG_succ::Chunk::build_graph_from_chunks(Gv.at(0)->get_k(),
                                                    { merged, });
}


std::vector<std::deque<TAlphabet>> get_last_added_nodes(const std::vector<const DBG_succ*> &Gv,
                                                        const std::vector<uint64_t> &kv) {
    std::vector<std::deque<TAlphabet>> last_added_nodes(DBG_succ::alph_size);
    // init last added nodes, if not starting from the beginning
    for (size_t i = 0; i < Gv.size(); i++) {
        // check whether we can merge the given graphs
        assert(Gv.at(i)->get_k() == Gv.at(0)->get_k()
                && "Graphs have different k-mer lengths - cannot be merged!\n");

        if (kv.at(i) < 2)
            continue;

        for (TAlphabet a = 0; a < DBG_succ::alph_size; a++) {
            uint64_t pred_pos = std::max(
                Gv.at(i)->pred_W(kv.at(i) - 1, a),
                Gv.at(i)->pred_W(kv.at(i) - 1, a + DBG_succ::alph_size)
            );
            if (pred_pos == 0)
                continue;

            auto curr_seq = Gv.at(i)->get_node_seq(pred_pos);

            if (!last_added_nodes[a].size()
                 || utils::colexicographically_greater(curr_seq, last_added_nodes[a]))
                last_added_nodes[a] = curr_seq;
        }
    }
    return last_added_nodes;
}

void merge_blocks(const std::vector<const DBG_succ*> &Gv,
                  std::vector<uint64_t> kv,
                  const std::vector<uint64_t> &nv,
                  DBG_succ::Chunk *chunk,
                  bool verbose) {

    assert(kv.size() == Gv.size());
    assert(nv.size() == Gv.size());

    auto last_added_nodes = get_last_added_nodes(Gv, kv);

    if (verbose) {
        std::cout << "Size of bins to merge: " << std::endl;
        for (size_t i = 0; i < Gv.size(); i++) {
            std::cout << nv.at(i) - kv.at(i) << std::endl;
        }
    }

    // Send parallel pointers running through each of the graphs. At each step, compare all
    // graph nodes at the respective positions with each other. Insert the lexicographically
    // smallest one into the common merge graph G (this).

    auto compare_edges = [](const std::deque<TAlphabet> &first,
                            const std::deque<TAlphabet> &second) {
        return utils::colexicographically_greater(second, first);
    };
    std::map<std::deque<TAlphabet>,
             std::vector<size_t>,
             decltype(compare_edges)> min_kmers(compare_edges);

    // find set of smallest pointers
    for (size_t i = 0; i < Gv.size(); ++i) {
        if (kv[i] < nv[i]) {
            auto seq = Gv[i]->get_node_seq(kv[i]);
            seq.push_front(Gv[i]->get_W(kv[i]) % DBG_succ::alph_size);
            min_kmers[seq].push_back(i);
        }
    }

    // keep track of how many nodes we added
    uint64_t added = 0;

    while (min_kmers.size()) {
        if (verbose && added > 0 && added % 10'000 == 0) {
            std::cout << "." << std::flush;
            if (added % 100'000 == 0) {
                std::cout << "added " << added;
                for (size_t i = 0; i < Gv.size(); i++)
                    std::cout << " - G" << i << ": edge " << kv.at(i)
                                             << "/" << Gv.at(i)->get_W().size();
                std::cout << std::endl;
            }
        }

        auto it = min_kmers.begin();

        auto seq1 = std::move(it->first);
        auto emptying_blocks = std::move(it->second);
        min_kmers.erase(it);

        TAlphabet val = seq1.front();
        seq1.pop_front();

        // check whether we already added a node whose outgoing edge points to the
        // same node as the current one
        TAlphabet next_in_W = val != DBG_succ::encode('$')
                                && utils::seq_equal(seq1, last_added_nodes[val], 1)
                              ? val + DBG_succ::alph_size
                              : val;

        bool remove_dummy_edge = false;

        // handle multiple outgoing edges
        if (chunk->size() > 0 && val != chunk->get_W_back() % DBG_succ::alph_size) {
            auto pred_node = last_added_nodes[chunk->get_W_back() % DBG_succ::alph_size];

            // compare the last two added nodes
            if (utils::seq_equal(seq1, pred_node)) {
                if (seq1.back() != DBG_succ::encode('$')
                        && chunk->get_W_back() == DBG_succ::encode('$')) {
                    remove_dummy_edge = true;
                    chunk->alter_W_back(next_in_W);
                } else {
                    chunk->alter_last_back(false);
                }
            }
        }
        if (!remove_dummy_edge)
            chunk->push_back(next_in_W, seq1.back(), true);

        last_added_nodes[val] = seq1;
        ++added;

        for (size_t i : emptying_blocks) {
            kv[i]++;
            if (kv[i] < nv[i]) {
                auto seq = Gv[i]->get_node_seq(kv[i]);
                seq.push_front(Gv[i]->get_W(kv[i]) % DBG_succ::alph_size);
                min_kmers[seq].push_back(i);
            }
        }
    }
}


} // namespace merge
