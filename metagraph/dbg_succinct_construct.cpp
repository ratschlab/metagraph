#include "dbg_succinct_construct.hpp"

#include <ips4o.hpp>

#include "kmer.hpp"
#include "dbg_succinct_chunk.hpp"
#include "utils.hpp"
#include "unix_tools.hpp"
#include "reads_filtering.hpp"


const size_t kMaxKmersChunkSize = 30'000'000;


template <class Vector>
void sort_and_remove_duplicates(Vector *kmers,
                                size_t num_threads,
                                size_t end_sorted = 0) {
    if (num_threads <= 3) {
        // sort
        ips4o::parallel::sort(kmers->data() + end_sorted,
                              kmers->data() + kmers->size(),
                              std::less<KMer>(), num_threads);
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
        ips4o::parallel::sort(kmers->data(), kmers->data() + kmers->size(),
                              std::less<KMer>(), num_threads);
    }
    // remove duplicates
    auto unique_end = std::unique(kmers->begin(), kmers->end());
    kmers->erase(unique_end, kmers->end());
}

template <class Vector>
void shrink_kmers(Vector *kmers,
                  size_t *end_sorted,
                  size_t num_threads,
                  bool verbose) {
    if (verbose) {
        std::cout << "Allocated capacity exceeded, filter out non-unique k-mers..."
                  << std::flush;
    }

    size_t prev_num_kmers = kmers->size();
    sort_and_remove_duplicates(kmers, num_threads, *end_sorted);
    *end_sorted = kmers->size();

    if (verbose) {
        std::cout << " done. Number of kmers reduced from " << prev_num_kmers
                                                  << " to " << *end_sorted << ", "
                  << (*end_sorted * sizeof(KMer) >> 20) << "Mb" << std::endl;
    }
}

template <class Array, class Vector>
void extend_kmer_storage(const Array &temp_storage,
                         Vector *kmers,
                         size_t *end_sorted,
                         size_t num_threads,
                         bool verbose,
                         std::mutex *mutex) {
    assert(mutex);

    // acquire the mutex to restrict the number of writing threads
    std::lock_guard<std::mutex> lock(*mutex);

    // shrink collected k-mers if the memory limit is exceeded
    if (kmers->size() + temp_storage.size() > kmers->capacity()) {
        shrink_kmers(kmers, end_sorted, num_threads, verbose);
        kmers->reserve(kmers->size()
                        + std::max(temp_storage.size(), kmers->size() / 2));
        if (kmers->size() + temp_storage.size() > kmers->capacity()) {
            std::cerr << "ERROR: Can't reallocate. Not enough memory" << std::endl;
        }
    }
    for (auto &kmer : temp_storage) {
        kmers->push_back(kmer);
    }
}

typedef std::function<void(const std::string&)> CallbackRead;

void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   KMerVector *kmers,
                   size_t *end_sorted,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex *mutex,
                   bool remove_redundant = true) {
    std::vector<KMer> temp_storage;
    temp_storage.reserve(1.1 * kMaxKmersChunkSize);

    generate_reads([&](const std::string &read) {
        utils::sequence_to_kmers(read, k, &temp_storage, suffix);

        if (temp_storage.size() < kMaxKmersChunkSize)
            return;

        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage, 1);
        }

        if (temp_storage.size() > 0.9 * kMaxKmersChunkSize) {
            extend_kmer_storage(temp_storage, kmers, end_sorted,
                                num_threads, verbose, mutex);
            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage, 1);
        }
        extend_kmer_storage(temp_storage, kmers, end_sorted,
                            num_threads, verbose, mutex);
    }
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
template <class Vector>
void recover_source_dummy_nodes(size_t k,
                                Vector *kmers,
                                size_t num_threads,
                                bool verbose) {
    // remove redundant dummy kmers inplace
    size_t cur_pos = 0;
    size_t end_sorted = kmers->size();

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
                                && KMer::compare_suffix(kmer, kmers->at(j), 1); ++j) {
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

        if (kmers->size() + k > kmers->capacity()) {
            if (verbose) {
                std::cout << "Allocated capacity exceeded,"
                          << " filter out non-unique k-mers..." << std::flush;
            }

            ips4o::parallel::sort(kmers->data() + end_sorted,
                                  kmers->data() + kmers->size(),
                                  std::less<KMer>(),
                                  num_threads);
            kmers->erase(
                std::unique(kmers->begin() + end_sorted, kmers->end()),
                kmers->end()
            );

            if (verbose) {
                std::cout << " done. Number of kmers: " << kmers->size() << ", "
                          << (kmers->size() * sizeof(KMer) >> 20) << "Mb" << std::endl;
            }
        }

        // anchor it to the dummy source node
        auto anchor_kmer = KMer::pack_kmer(std::vector<TAlphabet>(k + 1, 0), k + 1);
        for (size_t c = 2; c < k + 1; ++c) {
            KMer::update_kmer(k, kmers->at(i)[c], kmers->at(i)[c - 1], &anchor_kmer);

            kmers->emplace_back(anchor_kmer);
        }
    }
    std::copy(kmers->begin() + end_sorted, kmers->end(),
              kmers->begin() + cur_pos);
    kmers->resize(kmers->size() - end_sorted + cur_pos);

    sort_and_remove_duplicates(kmers, num_threads, cur_pos);
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
                                            double memory_preallocated,
                                            bool verbose)
      : k_(k),
        end_sorted_(0),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1,
                     std::max(static_cast<size_t>(1), num_threads_)),
        stored_reads_size_(0),
        verbose_(verbose) {
    assert(num_threads_ > 0);

    filter_suffix_encoded_.resize(filter_suffix.size());
    std::transform(filter_suffix.begin(), filter_suffix.end(),
                   filter_suffix_encoded_.begin(), DBG_succ::encode);

    kmers_.reserve(memory_preallocated / sizeof(KMer));

    if (filter_suffix == std::string(filter_suffix.size(), '$')) {
        kmers_.emplace_back(std::vector<TAlphabet>(k + 1, 0), k + 1);
    }
}

void KMerDBGSuccChunkConstructor::add_read(const std::string &sequence) {
    if (sequence.size() < k_)
        return;

    // put read into temporary storage
    stored_reads_size_ += sequence.size();
    reads_storage_.emplace_back(sequence);

    if (stored_reads_size_ < kMaxKmersChunkSize)
        return;

    // extract all k-mers from sequences accumulated in the temporary storage
    release_task_to_pool();

    assert(!stored_reads_size_);
    assert(!reads_storage_.size());
}

void KMerDBGSuccChunkConstructor::release_task_to_pool() {
    auto *current_reads_storage = new std::vector<std::string>();
    current_reads_storage->swap(reads_storage_);

    thread_pool_.enqueue(extract_kmers,
                         [current_reads_storage](CallbackRead callback) {
                             for (auto &&read : *current_reads_storage) {
                                 callback(std::move(read));
                             }
                             delete current_reads_storage;
                         },
                         k_, &kmers_, &end_sorted_, filter_suffix_encoded_,
                         num_threads_, verbose_, &mutex_, true);
    stored_reads_size_ = 0;
}

DBG_succ::Chunk* KMerDBGSuccChunkConstructor::build_chunk() {
    release_task_to_pool();
    thread_pool_.join();

    if (verbose_) {
        std::cout << "Reading data has finished" << std::endl;
        get_RAM();
        std::cout << "Sorting kmers and appending succinct"
                  << " representation from current bin...\t" << std::flush;
    }
    Timer timer;

    sort_and_remove_duplicates(&kmers_, num_threads_, end_sorted_);

    if (verbose_)
        std::cout << timer.elapsed() << "sec" << std::endl;

    if (!filter_suffix_encoded_.size()) {
        if (verbose_) {
            std::cout << "Reconstructing all required dummy source k-mers...\t"
                      << std::flush;
        }
        timer.reset();

        recover_source_dummy_nodes(k_, &kmers_, num_threads_, verbose_);

        if (verbose_)
            std::cout << timer.elapsed() << "sec" << std::endl;
    }

    DBG_succ::Chunk *result
        = DBG_succ::Chunk::build_from_kmers(k_, kmers_.data(), kmers_.size());

    kmers_.clear();

    return result;
}

void
KMerDBGSuccChunkConstructor
::add_reads(std::function<void(CallbackRead)> generate_reads) {
    thread_pool_.enqueue(extract_kmers, generate_reads,
                         k_, &kmers_, &end_sorted_,
                         filter_suffix_encoded_,
                         num_threads_, verbose_, &mutex_, true);
}


SuffixArrayDBGSuccConstructor::SuffixArrayDBGSuccConstructor(size_t k)
      : k_(k), data_("$") {}

void SuffixArrayDBGSuccConstructor::add_read(const std::string &read) {
    data_.append(read);
    data_.append("$");
}

// Implement SA construction and extract the kmers from the result
void SuffixArrayDBGSuccConstructor::build_graph(DBG_succ *graph) {
    DBG_succ::Chunk result;
    result.initialize_graph(graph);
}
