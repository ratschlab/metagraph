#include "dbg_succinct_merge.hpp"

#include <unistd.h>

#include "utils.hpp"
#include "dbg_succinct.hpp"
#include "config.hpp"


namespace merge {


pthread_mutex_t mutex_merge_result = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_bin_idx = PTHREAD_MUTEX_INITIALIZER;
// pthread_mutex_t mutex_annotate = PTHREAD_MUTEX_INITIALIZER;
pthread_attr_t attr;


/*
 * Helper function to determine the bin boundaries, given
 * a number of bins.
 */
std::vector<std::pair<uint64_t, uint64_t>> get_bins(const DBG_succ *G, uint64_t bins) {

    uint64_t nodes = G->rank_last(G->get_W().size() - 1);
    uint64_t orig_bins = bins;
    std::cerr << "working with " << orig_bins << " orig bins; "
                                     << nodes << " nodes" <<  std::endl;
    if (bins > nodes) {
        std::cerr << "[WARNING] There are max "
                  << nodes << " slots available for binning. Your current choice is "
                  << bins << " which will create "
                  << bins - nodes << " empty slots." << std::endl;
        bins = nodes;
    }

    std::vector<std::pair<uint64_t, uint64_t>> result;
    uint64_t binsize = (nodes + bins - 1) / bins;
    uint64_t thresh = (nodes - (bins * (nodes / bins))) * binsize;
    uint64_t pos = 1;
    for (uint64_t i = 0; i < nodes;) {
        if (i >= thresh) {
            binsize = nodes / bins;
        }
        result.push_back(
            std::make_pair(pos, G->select_last(std::min(nodes, i + binsize)))
        );
        pos = result.back().second + 1;
        i += binsize;
    }

    for (uint64_t i = bins; i < orig_bins; i++) {
        //result.push_back(std::make_pair(pos, pos));
        result.push_back(std::make_pair(1, 0));
    }

    std::cerr << "created " << result.size() << " bins" << std::endl;
    return result;
}


std::vector<std::pair<uint64_t, uint64_t>> get_bins_relative(
                                                const DBG_succ *G_from,
                                                const DBG_succ *G_to,
                                                std::vector<std::pair<uint64_t, uint64_t>> ref_bins,
                                                uint64_t first_pos,
                                                uint64_t last_pos) {

    std::vector<std::pair<uint64_t, uint64_t>> result;
    uint64_t pos = (first_pos == 0) ? 1 : G_from->colex_upper_bound(G_to->get_node_seq(first_pos)) + 1;
    uint64_t upper;
    for (size_t i = 0; i < ref_bins.size(); i++) {
        if (ref_bins.at(i).second == 0) { // this happens if we have more bins than nodes
            result.push_back(std::make_pair(0, 0));
        } else {
            upper = G_from->colex_upper_bound(G_to->get_node_seq(ref_bins.at(i).second));
            std::cerr << "ref bin " << ref_bins.at(i).second << " rel upper " << upper << std::endl;
            result.push_back(std::make_pair(pos, upper));
            pos = upper + 1;
        }
    }
    result.back().second = (last_pos == 0) ? G_from->get_W().size() - 1 : result.back().second;
    return result;
}


struct ParallelMergeContainer {
    std::vector<std::pair<uint64_t, uint64_t>> ref_bins;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> bins;
    std::vector<const DBG_succ*> result;
    std::vector<const DBG_succ*> graphs;
    unsigned int idx;
    unsigned int k;
    unsigned int bins_done;
    unsigned int first = 0;
    unsigned int last = 0;

    /* Helper function to subset the bins to the chunk
     * computed in the current distributed compute.
     */
    void subset_bins(unsigned int idx, unsigned int total, unsigned int bins_per_part){

        //std::cerr << "ref bins " << ref_bins.size() << " total " << total << " per part " << bins_per_part << std::endl;
        assert(ref_bins.size() == (total * bins_per_part));

        std::vector< std::pair<uint64_t, uint64_t> > new_ref_bins;
        //std::cerr << "min: " << binsize_min << " max: " << binsize_max
        //          << " thresh: " << threshold << " total: " << total << std::endl;

        /*size_t start, end;
        if (idx < threshold) {
            start = binsize_max * idx;
            end = (idx == (total - 1)) ? ref_bins.size() : binsize_max * (idx + 1);
        } else {
            start = (threshold * binsize_max) + ((idx - threshold) * binsize_min);
            end = (idx == (total - 1)) ? ref_bins.size() : (threshold * binsize_max) + ((idx - threshold + 1) * binsize_min);
        }*/
        size_t start = idx * bins_per_part;
        size_t end = (idx + 1) * bins_per_part;

        if (start > 0)
            first = ref_bins.at(start - 1).second;
        if (end < ref_bins.size())
            last = ref_bins.at(end).second;

        for (size_t i = start; i < end; i++) {
            new_ref_bins.push_back(ref_bins.at(i));
        }

        ref_bins = new_ref_bins;
    }

    void print_bins() {
        for (size_t ii = 0; ii < bins.size(); ii++) {
            std::cerr << "graph " << ii + 1 << std::endl;
            for (size_t i = 0; i < bins.at(ii).size(); i++)
                std::cerr << bins.at(ii).at(i).first << " - " << bins.at(ii).at(i).second << std::endl;
        }
    }

    /* Show an overview of the distribution of merging bin
     * sizes.
     */
    void get_bin_stats() {
        size_t min_bin = 0, max_bin = 0, total_bin = 0;

        size_t cum_size;
        for (size_t i = 0; i < ref_bins.size(); ++i) {
            cum_size = 0;
            for (size_t ii = 0; ii < bins.size(); ii++) {
                cum_size += (bins.at(ii).at(i).first == 0) ? 0 : bins.at(ii).at(i).second - bins.at(ii).at(i).first + 1;
            }
            if (cum_size > 0) {
                min_bin = (min_bin == 0) ? cum_size : std::min(min_bin, cum_size);
                max_bin = (max_bin == 0) ? cum_size : std::max(max_bin, cum_size);
            }
            total_bin += cum_size;
        }

        std::cout << std::endl;
        std::cout << "Total number of bins: " << ref_bins.size() << std::endl;
        std::cout << "Total size: " << total_bin << std::endl;
        std::cout << "Smallest bin: " << min_bin << std::endl;
        std::cout << "Largest bin: " << max_bin << std::endl;
        std::cout << "Average bin size: " << total_bin / ref_bins.size() << std::endl << std::endl;
    }
};


ParallelMergeContainer *merge_data;
// ParallelAnnotateContainer *anno_data = new ParallelAnnotateContainer();


/*
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used.
 */
void* parallel_merge_wrapper(void *config_) {
    Config *config = static_cast<Config *>(config_);

    unsigned int curr_idx;
    DBG_succ* graph;

    while (true) {
        pthread_mutex_lock (&mutex_bin_idx);
        if (merge_data->idx == merge_data->ref_bins.size()) {
            pthread_mutex_unlock (&mutex_bin_idx);
            break;
        } else {
            curr_idx = merge_data->idx;
            merge_data->idx++;
            pthread_mutex_unlock (&mutex_bin_idx);

            std::vector<uint64_t> kv;
            std::vector<uint64_t> nv;
            for (size_t i = 0; i < merge_data->graphs.size(); i++) {
                kv.push_back(merge_data->bins.at(i).at(curr_idx).first);
                nv.push_back(merge_data->bins.at(i).at(curr_idx).second + 1);
            }
            graph = merge::merge(merge_data->graphs, kv, nv);
            if (!graph) {
                std::cerr << "ERROR: Merge gone wrong" << std::endl;
                exit(1);
            }

            pthread_mutex_lock (&mutex_merge_result);
            merge_data->result.at(curr_idx) = graph;
            merge_data->bins_done++;
            if (config->verbose)
                std::cout << "finished bin " << curr_idx + 1
                          << " (" << merge_data->bins_done
                                  << "/" << merge_data->ref_bins.size()
                          << ")" << std::endl;
            pthread_mutex_unlock (&mutex_merge_result);
        }
    }
    pthread_exit((void*) 0);
}


DBG_succ* build_chunk(const std::vector<const DBG_succ*> &graphs, Config *config) {
    pthread_t *threads = NULL;
    merge_data = new ParallelMergeContainer();

    // get bins in graphs according to required threads
    if (config->verbose)
        std::cout << "Collecting reference bins" << std::endl;
    std::cerr << "parallel " << config->parallel
              << " per thread " << config->bins_per_thread
              << " parts total " << config->parts_total << std::endl;
    merge_data->ref_bins = get_bins(graphs.front(),
                                    config->parallel
                                        * config->bins_per_thread
                                        * config->parts_total);

    // only work on subset of the bins when requested
    if (config->parts_total > 1) {
        merge_data->subset_bins(config->part_idx,
                                config->parts_total,
                                config->parallel * config->bins_per_thread);
    }
    merge_data->bins.push_back(merge_data->ref_bins);

    if (config->verbose)
        std::cout << "Collecting relative bins" << std::endl;
    for (size_t i = 1; i < graphs.size(); i++) {
        std::cerr << "getting bins for " << i << ": " << config->fname[i] << std::endl;
        merge_data->bins.push_back(get_bins_relative(graphs.at(i), graphs.front(), merge_data->ref_bins,
                                                     merge_data->first, merge_data->last));
    }
    for (size_t i = 0; i < graphs.size(); i++) {
        for (size_t ii = 0; ii < merge_data->bins.at(i).size(); ii++) {
            if (merge_data->bins.at(i).at(ii).first > merge_data->bins.at(i).at(ii).second) {
               merge_data->bins.at(i).at(ii) = std::make_pair(graphs.at(i)->get_W().size(),
                                                              graphs.at(i)->get_W().size() - 1);
            }
        }
    }

    // print bin stats
    if (config->verbose) {
        merge_data->get_bin_stats();
    }

    // prepare data shared by threads
    merge_data->idx = 0;
    merge_data->k = graphs.front()->get_k();
    merge_data->graphs = graphs;
    for (size_t i = 0; i < merge_data->ref_bins.size(); i++)
        merge_data->result.push_back(NULL);
    merge_data->bins_done = 0;

    // create threads
    threads = new pthread_t[config->parallel];
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // do the work
    for (size_t tid = 0; tid < config->parallel; tid++) {
       pthread_create(&threads[tid], &attr, parallel_merge_wrapper, (void *) config);
       std::cerr << "starting thread " << tid << std::endl;
    }

    // join threads
    if (config->verbose)
        std::cout << "Waiting for threads to join" << std::endl;
    for (size_t tid = 0; tid < config->parallel; tid++) {
        pthread_join(threads[tid], NULL);
    }
    delete[] threads;

    // collect results
    std::cerr << "Collecting results" << std::endl;

    DBG_succ *graph = new DBG_succ(graphs.front()->get_k(), false);

    for (size_t i = 0; i < merge_data->result.size(); ++i) {
        if (merge_data->result.at(i)) {

            graph->verbose_cout("    adding ", merge_data->result.at(i)->W->size(), " edges\n");

            // handle last and W
            for (size_t j = 1, curr_pos = graph->W->size();
                                    j < merge_data->result.at(i)->W->size();
                                    ++j, ++curr_pos) {
                graph->W->insert(curr_pos, merge_data->result.at(i)->get_W(j));
                graph->last->insertBit(curr_pos, merge_data->result.at(i)->get_last(j));
            }

            graph->verbose_cout("new total edges: ", graph->W->size(), "\n");

            // handle F
            assert(graph->F.size() == merge_data->result.at(i)->F.size());
            for (size_t j = 0; j < graph->F.size(); ++j) {
                graph->F.at(j) += merge_data->result.at(i)->F.at(j);
            }

            delete merge_data->result.at(i);
        }
    }
    merge_data->result.clear();
    merge_data->bins_done = 0;

    delete merge_data;
    return graph;
}


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


DBG_succ* merge(const std::vector<const DBG_succ*> &Gv,
                std::vector<uint64_t> kv,
                std::vector<uint64_t> nv) {

    // Preliminarities
    DBG_succ *Gt = new DBG_succ(Gv.at(0)->get_k(), false);

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
        std::deque<TAlphabet> seq1 = Gv.at(curr_k)->get_node_seq(kv.at(curr_k));
        uint64_t val = Gv.at(curr_k)->get_W(kv.at(curr_k)) % Gt->alph_size;
        //std::cerr << "val: " << val << std::endl;
        //usleep(1000000);

        assert(val == smallest.second % Gt->alph_size);

        //std::cerr << "inserting into W" << std::endl;
        // check whether we already added a node whose outgoing edge points to the
        // same node as the current one
        auto it = last_added_nodes.find(smallest.second % Gt->alph_size);
        if (it != last_added_nodes.end() && utils::seq_equal(seq1, it->second, 1) && val != DBG_succ::encode('$')) {
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
                if (Gt->W->size() - 2 > 1 && Gt->get_W(Gt->W->size() - 2) == DBG_succ::encode('$')) {
                    Gt->update_F(Gv.at(curr_k)->get_node_last_char(kv.at(curr_k)), -1);
                    Gt->last->deleteBit(Gt->W->size() - 2);
                    Gt->W->remove(Gt->W->size() - 2);
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
