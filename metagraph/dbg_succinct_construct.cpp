#include "dbg_succinct_construct.hpp"

#include <parallel/algorithm>

#include "kmer.hpp"
#include "dbg_succinct_chunk.hpp"
#include "utils.hpp"


const size_t kNumBasepairsInTask = 1'000'000;

const size_t kMaxCounterSize = 10'000'000;


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
void sequence_to_kmers(std::vector<TAlphabet>&& seq,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix) {
    assert(k);
    assert(suffix.size() <= k);

    if (seq.size() < k + 1)
        return;

    // based on performance comparison
    // for KMer::pack_kmer and KMer::update_kmer
    if (suffix.size() > 1) {
        for (size_t i = 0; i < seq.size() - k; ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k] - suffix.size())) {
                kmers->emplace_back(&seq[i], k + 1);
            }
        }
    } else {
        // initialize and add the first kmer from sequence
        auto kmer = KMer::pack_kmer(seq.data(), k + 1);

        if (std::equal(suffix.begin(), suffix.end(),
                       &seq[k] - suffix.size())) {
            kmers->emplace_back(kmer);
        }

        // add all other kmers
        for (size_t i = 1; i < seq.size() - k; ++i) {
            KMer::update_kmer(k, seq[i + k], seq[i + k - 1], &kmer);

            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k] - suffix.size())) {
                kmers->emplace_back(kmer);
            }
        }
    }
}

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
void sequence_to_kmers(const std::string &sequence,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix) {
    assert(k);
    assert(suffix.size() <= k);

    if (sequence.size() < k)
        return;

    // encode sequence
    size_t dummy_prefix_size = suffix.size() > 0 ? k : 1;

    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1);

    for (size_t i = 0; i < dummy_prefix_size; ++i) {
        seq[i] = DBG_succ::encode('$');
    }
    std::transform(sequence.begin(), sequence.end(),
                   &seq[dummy_prefix_size], DBG_succ::encode);
    seq.back() = DBG_succ::encode('$');

    sequence_to_kmers(std::move(seq), k, kmers, suffix);
}

void sort_and_remove_duplicates(std::vector<KMer> *kmers,
                                size_t num_threads,
                                size_t end_sorted = 0);

void shrink_kmers(std::vector<KMer> *kmers,
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

// takes the ownership of the allocated sequence and releases when finishes
void extend_kmer_storage(const std::vector<KMer> &temp_storage,
                         std::vector<KMer> *kmers,
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
    // try {
    for (auto &kmer : temp_storage) {
        kmers->push_back(kmer);
    }
    // } catch (...) {
    //     std::cerr << "ERROR: Not enough memory."
    //               << " Try to increase the memory limit." << std::endl;
    //     exit(1);
    // }
}

// takes the ownership of the allocated sequence and releases when finishes
void sequence_to_kmers_parallel(std::vector<std::string> *reads,
                                size_t k,
                                std::vector<KMer> *kmers,
                                size_t *end_sorted,
                                const std::vector<TAlphabet> &suffix,
                                size_t num_threads,
                                bool verbose,
                                std::mutex *mutex,
                                bool remove_redundant) {
    assert(mutex);

    // parallel mode
    std::vector<KMer> temp_storage;

    size_t num_basepairs = 0;
    for (auto &read : *reads) {
        num_basepairs += read.size() + 2;
    }
    temp_storage.reserve(num_basepairs);

    for (const auto &read : *reads) {
        sequence_to_kmers(read, k, &temp_storage, suffix);
    }
    delete reads;

    if (remove_redundant) {
        sort_and_remove_duplicates(&temp_storage, 1, 0);
    }

    extend_kmer_storage(temp_storage, kmers, end_sorted,
                        num_threads, verbose, mutex);
}

void sort_and_remove_duplicates(std::vector<KMer> *kmers,
                                size_t num_threads,
                                size_t end_sorted) {
    if (num_threads <= 3) {
        // sort
        omp_set_num_threads(num_threads);
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
        omp_set_num_threads(num_threads);
        __gnu_parallel::sort(kmers->data(), kmers->data() + kmers->size());
    }
    // remove duplicates
    auto unique_end = std::unique(kmers->begin(), kmers->end());
    kmers->erase(unique_end, kmers->end());
}

void recover_source_dummy_nodes(size_t k,
                                std::vector<KMer> *kmers,
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

        if (kmers->size() + k > kmers->capacity()) {
            if (verbose) {
                std::cout << "Allocated capacity exceeded,"
                          << " filter out non-unique k-mers..." << std::flush;
            }

            omp_set_num_threads(num_threads);
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
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1),
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

    if (stored_reads_size_ < kNumBasepairsInTask)
        return;

    // extract all k-mers from sequences accumulated in the temporary storage
    release_task_to_pool();

    assert(!stored_reads_size_);
    assert(!reads_storage_.size());
}

void KMerDBGSuccChunkConstructor::release_task_to_pool() {
    auto *current_reads_storage = new std::vector<std::string>();
    current_reads_storage->swap(reads_storage_);
    stored_reads_size_ = 0;

    thread_pool_.enqueue(sequence_to_kmers_parallel, current_reads_storage,
                         k_, &kmers_, &end_sorted_, filter_suffix_encoded_,
                         num_threads_, verbose_, &mutex_, true);
}

DBG_succ::Chunk* KMerDBGSuccChunkConstructor::build_chunk() {
    release_task_to_pool();
    thread_pool_.join();
    sort_and_remove_duplicates(&kmers_, num_threads_, end_sorted_);

    if (!filter_suffix_encoded_.size()) {
        recover_source_dummy_nodes(k_, &kmers_, num_threads_, verbose_);
    }

    DBG_succ::Chunk *result = DBG_succ::VectorChunk::build_from_kmers(k_, &kmers_);
    kmers_.clear();

    return result;
}

typedef std::function<void(const std::string&)> CallbackRead;

void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   std::vector<KMer> *kmers,
                   size_t *end_sorted,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex *mutex) {
    std::vector<KMer> temp_storage;

    generate_reads([&](const std::string &read) {
        sequence_to_kmers(read, k, &temp_storage, suffix);

        if (temp_storage.size() < 20'000'000)
            return;

        sort_and_remove_duplicates(&temp_storage, 1, 0);

        if (temp_storage.size() < 18'000'000)
            return;

        extend_kmer_storage(temp_storage, kmers, end_sorted,
                            num_threads, verbose, mutex);
        temp_storage.clear();
    });

    sort_and_remove_duplicates(&temp_storage, 1, 0);
    extend_kmer_storage(temp_storage, kmers, end_sorted,
                        num_threads, verbose, mutex);
}

struct KMerHash {
    size_t operator()(const KMer &kmer) const {
        const uint32_t *ker_ptr = reinterpret_cast<const uint32_t*>(&kmer);
        // computes the hash of a k-mer using a variant
        // of the Fowler-Noll-Vo hash function
        size_t result = ker_ptr[0];

        for (size_t i = 1; i < sizeof(KMer) / sizeof(uint32_t); ++i) {
            result = (result * 16777619) ^ ker_ptr[i];
        }

        return result;
    }
};

// use a custom allocator since the default one is too slow
typedef std::unordered_map<KMer, uint32_t, KMerHash, std::equal_to<KMer>,
                           utils::plalloc<std::pair<const KMer, uint32_t>>> Counter;

void move_kmers_to_storage(Counter *counter,
                           std::vector<KMer> *kmers,
                           size_t *end_sorted,
                           size_t noise_kmer_frequency,
                           size_t num_threads,
                           bool verbose,
                           std::mutex *mutex) {
    std::vector<KMer> temp_storage;

    size_t counter_sum = 0;

    for (auto it = counter->begin(); it != counter->end(); ++it) {
        counter_sum += it->second;
        if (it->second > noise_kmer_frequency)
            temp_storage.push_back(it->first);
    }

    if (verbose) {
        std::cout << "\nFiltering out the k-mers collected...\n";
        std::cout << "Total k-mers:    " << counter_sum << std::endl;
        std::cout << "Distinct k-mers: " << counter->size() << std::endl;
        std::cout << "Filtered k-mers: " << temp_storage.size() << std::endl;
    }
    counter->clear();

    extend_kmer_storage(temp_storage, kmers, end_sorted,
                        num_threads, verbose, mutex);
}

void count_kmers(std::function<void(CallbackRead)> generate_reads,
                 size_t k,
                 std::vector<KMer> *kmers,
                 size_t *end_sorted,
                 const std::vector<TAlphabet> &suffix,
                 size_t noise_kmer_frequency,
                 size_t num_threads,
                 bool verbose,
                 std::mutex *mutex) {
    if (noise_kmer_frequency == 0) {
        extract_kmers(generate_reads, k, kmers, end_sorted, suffix,
                      num_threads, verbose, mutex);
        return;
    }

    Counter counter;
    counter.reserve(kMaxCounterSize);
    std::vector<KMer> temp_storage;

    generate_reads([&](const std::string &read) {
        sequence_to_kmers(read, k, &temp_storage, suffix);

        for (const auto &kmer : temp_storage) {
            counter[kmer]++;
        }
        temp_storage.resize(0);

        if (counter.size() > kMaxCounterSize) {
            move_kmers_to_storage(&counter, kmers, end_sorted,
                                  noise_kmer_frequency,
                                  num_threads, verbose, mutex);
            counter.reserve(kMaxCounterSize);
        }
    });
    move_kmers_to_storage(&counter, kmers, end_sorted,
                          noise_kmer_frequency * counter.size() / kMaxCounterSize,
                          num_threads, verbose, mutex);
}

void KMerDBGSuccChunkConstructor::add_reads(std::function<void(CallbackRead)> generate_reads,
                                            size_t noise_kmer_frequency) {
    thread_pool_.enqueue(count_kmers, generate_reads,
                         k_, &kmers_, &end_sorted_,
                         filter_suffix_encoded_, noise_kmer_frequency,
                         num_threads_, verbose_, &mutex_);
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
