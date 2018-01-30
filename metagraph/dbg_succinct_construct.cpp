#include "dbg_succinct_construct.hpp"

#include <parallel/algorithm>

#include "kmer.hpp"
#include "dbg_succinct_chunk.hpp"
#include "utils.hpp"


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
    auto kmer = KMer::pack_kmer(seq.data(), k + 1);

    if (std::equal(suffix.begin(), suffix.end(),
                   seq.data() + k - suffix.size())) {
        kmers->emplace_back(kmer);
    }

    // add all other kmers
    for (size_t i = 1; i < seq.size() - k; ++i) {
        update_kmer(k, seq[i + k], seq[i + k - 1], &kmer);

        if (std::equal(suffix.begin(), suffix.end(),
                       seq.data() + i + k - suffix.size())) {
            kmers->emplace_back(kmer);
        }
    }
}

void sort_and_remove_duplicates(std::vector<KMer> *kmers, size_t end_sorted = 0) {
    if (omp_get_num_threads() <= 3) {
        // sort
        __gnu_parallel::sort(kmers->data() + end_sorted,
                             kmers->data() + kmers->size());
        kmers->erase(std::unique(kmers->begin() + end_sorted, kmers->end()),
                     kmers->end());

        // merge two sorted arrays
#if 0
        std::__merge_without_buffer(kmers->data(),
                                    kmers->data() + end_sorted,
                                    kmers->data() + kmers->size(),
                                    end_sorted, kmers->size() - end_sorted,
                                    __gnu_cxx::__ops::__iter_less_iter());
#else
        std::inplace_merge(kmers->data(),
                           kmers->data() + end_sorted,
                           kmers->data() + kmers->size());
#endif
    } else {
        // sort
        __gnu_parallel::sort(kmers->data(), kmers->data() + kmers->size());
    }
    // remove duplicates
    auto unique_end = std::unique(kmers->begin(), kmers->end());
    kmers->erase(unique_end, kmers->end());
}

void recover_source_dummy_nodes(size_t k,
                                std::vector<KMer> *kmers,
                                size_t max_num_kmers,
                                bool verbose) {
    // remove redundant dummy kmers inplace
    size_t cur_pos = 0;
    size_t end_sorted = kmers->size();

    kmers->emplace_back(KMer::pack_kmer(std::vector<TAlphabet>(k + 1, 0), k + 1));

    for (size_t i = 0; i < end_sorted; ++i) {
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
        for (size_t j = i + 1; j < end_sorted
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

        if (kmers->size() + k > max_num_kmers) {
            if (verbose) {
                std::cout << "Memory limit exceeded, filter out non-unique k-mers..." << std::flush;
            }

            __gnu_parallel::sort(kmers->data() + end_sorted,
                                 kmers->data() + kmers->size());
            kmers->erase(
                std::unique(kmers->begin() + end_sorted, kmers->end()),
                kmers->end()
            );

            if (verbose) {
                std::cout << " done. Number of kmers: " << kmers->size() << ", "
                          << (kmers->size() * sizeof(KMer) >> 20) << "Mb" << std::endl;
            }

            if (kmers->size() + k > max_num_kmers - max_num_kmers / 100) {
                std::cerr << "ERROR: Not enough memory."
                          << " Try to increase the memory limit." << std::endl;
                exit(1);
            }
        }

        // anchor it to the dummy source node
        auto anchor_kmer = KMer::pack_kmer(std::vector<TAlphabet>(k + 1, 0), k + 1);
        for (size_t c = 2; c < k + 1; ++c) {
            update_kmer(k, kmers->at(i)[c], kmers->at(i)[c - 1], &anchor_kmer);

            kmers->emplace_back(anchor_kmer);
        }
    }
    std::copy(kmers->begin() + end_sorted, kmers->end(),
              kmers->begin() + cur_pos);
    kmers->resize(kmers->size() - end_sorted + cur_pos);

    sort_and_remove_duplicates(kmers, cur_pos);
}


void KMerDBGSuccConstructor::build_graph(DBG_succ *graph) {
    // build the graph chunk from kmers
    auto chunk = constructor_.build_chunk();
    // initialize graph from the chunk built
    chunk->initialize_graph(graph);
    delete chunk;
}


KMerDBGSuccChunkConstructor::KMerDBGSuccChunkConstructor(
                                            size_t k,
                                            const std::string &filter_suffix,
                                            size_t num_threads,
                                            double memory_available,
                                            bool verbose)
      : k_(k),
        end_sorted_(0),
        num_threads_(num_threads),
        max_num_kmers_(memory_available / sizeof(KMer)),
        verbose_(verbose) {
    filter_suffix_encoded_.resize(filter_suffix.size());
    std::transform(filter_suffix.begin(), filter_suffix.end(),
                   filter_suffix_encoded_.begin(), DBG_succ::encode);
    if (max_num_kmers_ == 0) {
        max_num_kmers_ = static_cast<size_t>(-1);
    } else {
        kmers_.reserve(max_num_kmers_);
    }
    omp_set_num_threads(std::max(static_cast<int>(num_threads_), 1));
}

void KMerDBGSuccChunkConstructor::add_read(const std::string &sequence) {
    if (kmers_.size() + sequence.size() > max_num_kmers_) {
        if (verbose_) {
            std::cout << "Memory limit exceeded, filter out non-unique k-mers..." << std::flush;
        }

        sort_and_remove_duplicates(&kmers_, end_sorted_);
        end_sorted_ = kmers_.size();

        if (verbose_) {
            std::cout << " done. Number of kmers: " << end_sorted_ << ", "
                      << (end_sorted_ * sizeof(KMer) >> 20) << "Mb" << std::endl;
        }

        if (kmers_.size() + sequence.size() > max_num_kmers_ - max_num_kmers_ / 50) {
            std::cerr << "ERROR: Not enough memory."
                      << " Try to increase the memory limit." << std::endl;
            exit(1);
        }
    }

    // add all k-mers of seq to the graph
    sequence_to_kmers(sequence, k_, &kmers_, filter_suffix_encoded_);
}

DBG_succ::Chunk* KMerDBGSuccChunkConstructor::build_chunk() {
    sort_and_remove_duplicates(&kmers_, end_sorted_);

    if (!filter_suffix_encoded_.size()) {
        recover_source_dummy_nodes(k_, &kmers_, max_num_kmers_, verbose_);
    }

    DBG_succ::Chunk *result = DBG_succ::VectorChunk::build_from_kmers(k_, &kmers_);
    kmers_.clear();

    return result;
}


SuffixArrayDBGSuccConstructor::SuffixArrayDBGSuccConstructor(size_t k)
      : k_(k), data_("$") {}

void SuffixArrayDBGSuccConstructor::add_read(const std::string &read) {
    data_.append(read);
    data_.append("$");
}

// Implement SA construction and extract the kmers from the result
void SuffixArrayDBGSuccConstructor::build_graph(DBG_succ *graph) {
    DBG_succ::VectorChunk result;
    result.initialize_graph(graph);
}
