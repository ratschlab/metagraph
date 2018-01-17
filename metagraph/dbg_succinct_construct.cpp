#include "dbg_succinct_construct.hpp"

#include <parallel/algorithm>

#include "kmer.hpp"
#include "dbg_succinct_chunk.hpp"


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
void sequence_to_kmers(const std::string &sequence,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix = {}) {
    if (sequence.size() < k)
        return;

    // encode sequence
    size_t dummy_prefix_size = suffix.size() > 0 ? k + 1 : 1;
    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1);
    for (size_t i = 0; i < dummy_prefix_size; ++i) {
        seq[i] = DBG_succ::encode('$');
    }
    std::transform(sequence.begin(), sequence.end(),
                   &seq[dummy_prefix_size], DBG_succ::encode);
    seq.back() = DBG_succ::encode('$');

    // initialize and add the first kmer from sequence
    KMer kmer(seq.data(), k + 1);

    if (std::equal(suffix.begin(), suffix.end(),
                   seq.data() + k - suffix.size())) {
        kmers->emplace_back(kmer);
    }

    // add all other kmers
    for (size_t i = 1; i < seq.size() - k; ++i) {
        kmer.update(k, seq[i + k]);

        if (std::equal(suffix.begin(), suffix.end(),
                       seq.data() + i + k - suffix.size())) {
            kmers->emplace_back(kmer);
        }
    }
}

void sort_and_remove_duplicates(std::vector<KMer> *kmers) {
    // sort
    __gnu_parallel::sort(kmers->data(), kmers->data() + kmers->size());

    // remove duplicates
    auto unique_end = std::unique(kmers->begin(), kmers->end());
    kmers->erase(unique_end, kmers->end());
}

void recover_source_dummy_nodes(size_t k, std::vector<KMer> *kmers) {
    // remove redundant dummy kmers inplace
    size_t cur_pos = 0;
    std::vector<KMer> prev_dummy_kmers {
        KMer(std::vector<TAlphabet>(k + 1, 0), k + 1)
    };

    for (size_t i = 0; i < kmers->size(); ++i) {
        const KMer &kmer = kmers->at(i);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        TAlphabet edge_label;

        // check if it's not a source dummy kmer
        if (kmer[1] > 0 || (edge_label = kmer[0]) == 0) {
            kmers->at(cur_pos++) = kmer;
            continue;
        }

        bool redundant = false;
        for (size_t j = i + 1; j < kmers->size()
                                    && KMer::compare_kmer_suffix(kmer, kmers->at(j), 1); ++j) {
            if (edge_label == kmers->at(j)[0]) {
                // This source dummy kmer is redundant and has to be erased
                redundant = true;
                break;
            }
        }
        if (redundant)
            continue;

        // leave this dummy kmer in the list
        kmers->at(cur_pos++) = kmer;

        // anchor it to the dummy source node
        KMer anchor_kmer(std::vector<TAlphabet>(k + 1, 0), k + 1);
        for (size_t c = 2; c < k + 1; ++c) {
            anchor_kmer.update(k, kmer[c]);
            prev_dummy_kmers.emplace_back(anchor_kmer);
        }
    }
    kmers->resize(cur_pos);

    sort_and_remove_duplicates(&prev_dummy_kmers);

    std::vector<KMer> proper_kmers(kmers->size() + prev_dummy_kmers.size());
    kmers->swap(proper_kmers);
    __gnu_parallel::merge(proper_kmers.begin(), proper_kmers.end(),
                          prev_dummy_kmers.begin(), prev_dummy_kmers.end(),
                          kmers->data());
}


KMerDBGSuccConstructor::KMerDBGSuccConstructor(size_t k, size_t num_threads)
      : k_(k), num_threads_(num_threads) {}

void KMerDBGSuccConstructor::add_read(const std::string &read) {
    sequence_to_kmers(read, k_, &kmers_);
}

void KMerDBGSuccConstructor::add_reads(const std::vector<std::string> &reads) {
    // break the sequences down into kmers
    for (const auto &sequence : reads) {
        add_read(sequence);
    }
}

void KMerDBGSuccConstructor::build_graph(DBG_succ *graph) {
    assert(k_ == graph->get_k());

    omp_set_num_threads(std::max(static_cast<int>(num_threads_), 1));

    sort_and_remove_duplicates(&kmers_);

    recover_source_dummy_nodes(k_, &kmers_);

    // build the graph chunk from kmers
    auto chunk = DBG_succ::VectorChunk::build_from_kmers(k_, &kmers_);
    // initialize graph from the chunk built
    chunk->initialize_graph(graph);
    delete chunk;
}


KMerDBGSuccChunkConstructor::KMerDBGSuccChunkConstructor(
                                            size_t k,
                                            const std::string &filter_suffix,
                                            size_t num_threads
) : k_(k), num_threads_(num_threads) {
    filter_suffix_encoded_.resize(filter_suffix.size());
    std::transform(filter_suffix.begin(), filter_suffix.end(),
                   filter_suffix_encoded_.begin(), DBG_succ::encode);
}

void KMerDBGSuccChunkConstructor::add_read(const std::string &sequence) {
    // add all k-mers of seq to the graph
    // TODO: This reserve makes the program super slow... why?..
    // kmers.reserve(kmers.size() + read_stream->seq.l);
    sequence_to_kmers(sequence, k_, &kmers_, filter_suffix_encoded_);
}

DBG_succ::Chunk* KMerDBGSuccChunkConstructor::build_chunk() {
    omp_set_num_threads(std::max(static_cast<int>(num_threads_), 1));

    sort_and_remove_duplicates(&kmers_);

    if (!filter_suffix_encoded_.size())
        recover_source_dummy_nodes(k_, &kmers_);

    DBG_succ::Chunk *result = DBG_succ::VectorChunk::build_from_kmers(k_, &kmers_);
    kmers_.clear();

    return result;
}
