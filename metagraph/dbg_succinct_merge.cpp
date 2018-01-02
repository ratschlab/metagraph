#include "dbg_succinct_merge.hpp"

#include <unistd.h>
#include <thread>
#include <mutex>

#include "utils.hpp"
#include "dbg_succinct.hpp"
#include "config.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::NumberSerialisation;


namespace merge {

/**
 * Helper function to determine the bin boundaries, given
 * a number of bins basen on the total number of nodes in graph.
 * [v[0], v[1]), [v[1], v[2]), ...
 */
std::vector<uint64_t> get_bins(const DBG_succ &G, uint64_t num_bins) {
    assert(num_bins > 0);

    uint64_t num_nodes = G.rank_last(G.get_W().size() - 1);

    if (G.verbose && num_bins > num_nodes) {
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


/*
// TODO: Do we need this version? Probably more memory efficient but it's much slower.
class dynamic_graph_chunk : public graph_chunk {
  public:
    void push_back(TAlphabet W, TAlphabet F, bool last) {
        W_.insert(W_.size(), W);
        for (TAlphabet a = F + 1; a < F_.size(); ++a) {
            F_[a]++;
        }
        last_.insertBit(last_.size(), last);
    }

    TAlphabet get_W_back() const { return W_[W_.size() - 1]; }

    void alter_W_back(TAlphabet W) { W_.set(W_.size() - 1, W); }

    void alter_last_back(bool last) { last_.set(last_.size() - 1, last); }

    uint64_t size() const { return W_.size(); }

    void extend(const graph_chunk &other) {
        extend(dynamic_cast<const dynamic_graph_chunk&>(other));
    }

    void extend(const dynamic_graph_chunk &other) {
        for (uint64_t j = 0; j < other.W_.size(); ++j) {
            W_.insert(W_.size(), other.W_[j]);
            last_.insertBit(last_.size(), other.last_[j]);
        }

        assert(F_.size() == other.F_.size());
        for (size_t p = 0; p < other.F_.size(); ++p) {
            F_[p] += other.F_[p];
        }
    }

    void initialize_graph(DBG_succ *graph) {
        delete graph->W;
        graph->W = new wavelet_tree_dyn(4, W_);

        delete graph->last;
        graph->last = new bit_vector_dyn(last_);

        graph->F = F_;
    }

    bool load(const std::string &infbase) {
        F_.resize(0);

        try {
            std::ifstream instream(infbase + ".F.dbg");
            std::string cur_line;

            while (std::getline(instream, cur_line)) {
                F_.push_back(std::stoull(cur_line));
            }
            instream.close();

            if (F_.size() != DBG_succ::alph_size)
                return false;

            // load W and last arrays
            std::ifstream instream_W(infbase + ".W.dbg");
            std::ifstream instream_l(infbase + ".l.dbg");
            return W_.deserialise(instream_W) && last_.deserialise(instream_l);
        } catch (...) {
            return false;
        }
    }

    void serialize(const std::string &outbase) const {
        std::ofstream outstream(outbase + ".W.dbg");
        W->serialise(outstream);
        outstream.close();

        // write last array
        outstream.open(outbase + ".l.dbg");
        last->serialise(outstream);
        outstream.close();

        // write F values and k
        outstream.open(outbase + ".F.dbg");
        for (size_t i = 0; i < F.size(); ++i) {
            outstream << F.at(i) << "\n";
        }
        outstream.close();
    }

  private:
    wavelet_tree_dyn W_ = wavelet_tree_dyn(4);
    bit_vector_dyn last_;
    std::vector<uint64_t> F_ = std::vector<uint64_t>(DBG_succ::alph_size, 0);
};
*/

class vector_graph_chunk : public graph_chunk {
  public:
    void push_back(TAlphabet W, TAlphabet F, bool last) {
        W_.push_back(W);
        for (TAlphabet a = F + 1; a < F_.size(); ++a) {
            F_[a]++;
        }
        last_.push_back(last);
    }

    TAlphabet get_W_back() const { return W_.back(); }

    void alter_W_back(TAlphabet W) { W_.back() = W; }

    void alter_last_back(bool last) { last_.back() = last; }

    uint64_t size() const { return W_.size(); }

    void extend(const graph_chunk &other) {
        extend(dynamic_cast<const vector_graph_chunk&>(other));
    }

    void extend(const vector_graph_chunk &other) {
        W_.insert(W_.end(), other.W_.begin(), other.W_.end());
        last_.insert(last_.end(), other.last_.begin(), other.last_.end());

        assert(F_.size() == other.F_.size());
        for (size_t p = 0; p < other.F_.size(); ++p) {
            F_[p] += other.F_[p];
        }
    }

    void initialize_graph(DBG_succ *graph) {
        delete graph->W;
        graph->W = new wavelet_tree_dyn(4, W_);

        delete graph->last;
        graph->last = new bit_vector_dyn(last_);

        graph->F = F_;
    }

    bool load(const std::string &infbase) {
        try {
            std::ifstream instream_W(infbase + ".W.dbg");
            W_ = NumberSerialisation::deserialiseNumberVector<TAlphabet>(instream_W);
            instream_W.close();

            std::ifstream instream_l(infbase + ".l.dbg");
            last_ = NumberSerialisation::deserialiseNumberVector<bool>(instream_l);
            instream_l.close();

            std::ifstream instream_F(infbase + ".F.dbg");
            F_ = NumberSerialisation::deserialiseNumberVector<uint64_t>(instream_F);
            instream_F.close();

            return F_.size() == DBG_succ::alph_size;
        } catch (...) {
            return false;
        }
    }

    void serialize(const std::string &outbase) const {
        std::ofstream outstream(outbase + ".W.dbg");
        NumberSerialisation::serialiseNumberVector(outstream, W_);
        outstream.close();

        outstream.open(outbase + ".l.dbg");
        NumberSerialisation::serialiseNumberVector(outstream, last_);
        outstream.close();

        outstream.open(outbase + ".F.dbg");
        NumberSerialisation::serialiseNumberVector(outstream, F_);
        outstream.close();
    }

  private:
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_ = std::vector<uint64_t>(DBG_succ::alph_size, 0);
};


void merge_blocks(const std::vector<const DBG_succ*> &Gv,
                  std::vector<uint64_t> kv,
                  const std::vector<uint64_t> &nv,
                  graph_chunk *chunk);

/**
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used.
 */
void parallel_merge_wrapper(const std::vector<const DBG_succ*> &graphs,
                            const std::vector<std::vector<uint64_t>> &bins,
                            ParallelMergeContext *context,
                            std::vector<graph_chunk*> *chunks) {
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
        merge_blocks(graphs, kv, nv, chunks->at(curr_idx));
    }
}


graph_chunk* build_chunk(const std::vector<const DBG_succ*> &graphs,
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
        get_bins(*graphs.front(),
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

    std::vector<graph_chunk*> blocks(bins.front().size() - 1);
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = new vector_graph_chunk();
    }

    for (size_t tid = 0; tid < num_threads; tid++) {
        threads.emplace_back(parallel_merge_wrapper, graphs,
                                                     bins,
                                                     &context,
                                                     &blocks);
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

    graph_chunk *chunk = new vector_graph_chunk();

    for (uint64_t i = 0; i < blocks.size(); ++i) {
        chunk->extend(*blocks[i]);
        delete blocks[i];
    }

    return chunk;
}


/**
 * Merge graph chunks from the vector passed and release the chunks afterwards.
 * If the null pointers are passed,
 * load chunks from files "filenamebase.<chunk_idx>_<num_chunks>".
 */
DBG_succ* merge_chunks(size_t k,
                       const std::vector<graph_chunk*> &graph_chunks,
                       const std::string &filenamebase) {

    DBG_succ *graph = new DBG_succ(k);
    if (!graph_chunks.size())
        return graph;

    graph_chunk *merged = new vector_graph_chunk();
    merged->push_back(0, DBG_succ::alph_size, 0);

    for (uint64_t f = 0; f < graph_chunks.size(); f++) {
        graph_chunk *chunk = graph_chunks[f];

        // load from file if undefined
        if (!chunk) {
            std::string filename = filenamebase
                                    + "." + std::to_string(f)
                                    + "_" + std::to_string(graph_chunks.size());

            chunk = new vector_graph_chunk();
            if (!chunk->load(filename)) {
                std::cerr << "ERROR: input file "
                          << filename << " corrupted" << std::endl;
                delete graph;
                delete merged;
                delete chunk;
                return NULL;
            }
        }
        merged->extend(*chunk);
        delete chunk;
    }

    merged->initialize_graph(graph);
    delete merged;

    return graph;
}


DBG_succ* merge(const std::vector<const DBG_succ*> &Gv) {
    std::vector<uint64_t> kv;
    std::vector<uint64_t> nv;

    for (size_t i = 0; i < Gv.size(); ++i) {
        kv.push_back(1);
        nv.push_back(Gv[i]->get_W().size());
    }

    vector_graph_chunk merged;
    merged.push_back(0, DBG_succ::alph_size, 0);

    merge_blocks(Gv, kv, nv, &merged);

    DBG_succ *graph = new DBG_succ(Gv.at(0)->get_k());
    merged.initialize_graph(graph);

    return graph;
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
                  graph_chunk *chunk) {

    assert(kv.size() == Gv.size());
    assert(nv.size() == Gv.size());

    auto verbose = Gv.at(0)->verbose;

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

    // keep track of how many nodes we added
    uint64_t added = 0;

    while (true) {
        if (verbose && added > 0 && added % 1000 == 0) {
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
        auto smallest = utils::smallest_nodes(Gv, kv, nv);

        auto it = std::find(smallest.begin(), smallest.end(), true);
        if (it == smallest.end())
            break;

        size_t i = it - smallest.begin();

        auto seq1 = Gv.at(i)->get_node_seq(kv.at(i));

        TAlphabet val = Gv.at(i)->get_W(kv.at(i)) % DBG_succ::alph_size;

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

        uint64_t updated = 0;
        for (size_t i = 0; i < Gv.size(); i++) {
            if (smallest.at(i)) {
                updated++;
                kv.at(i)++;
            }
        }
        if (updated == 0)
            break;
    }
}


} // namespace merge
