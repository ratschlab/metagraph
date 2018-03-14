#include "dbg_succinct_construct.hpp"

#include <parallel/algorithm>
#include <hopscotch_map.h>

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

template <class Array>
void extend_kmer_storage(const Array &temp_storage,
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

typedef std::function<void(const std::string&)> CallbackRead;

void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   std::vector<KMer> *kmers,
                   size_t *end_sorted,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex *mutex,
                   bool remove_redundant = true,
                   size_t num_appended_kmers = kMaxCounterSize) {
    std::vector<KMer> temp_storage;
    temp_storage.reserve(1.1 * num_appended_kmers);

    generate_reads([&](const std::string &read) {
        sequence_to_kmers(read, k, &temp_storage, suffix);

        if (temp_storage.size() < num_appended_kmers)
            return;

        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage, num_threads);
        }

        if (temp_storage.size() > 0.9 * num_appended_kmers) {
            extend_kmer_storage(temp_storage, kmers, end_sorted,
                                num_threads, verbose, mutex);
            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage, num_threads);
        }
        extend_kmer_storage(temp_storage, kmers, end_sorted,
                            num_threads, verbose, mutex);
    }
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

    thread_pool_.enqueue(extract_kmers,
                         [current_reads_storage](CallbackRead callback) {
                             for (auto &&read : *current_reads_storage) {
                                 callback(std::move(read));
                             }
                             delete current_reads_storage;
                         },
                         k_, &kmers_, &end_sorted_, filter_suffix_encoded_,
                         num_threads_, verbose_, &mutex_, true,
                         std::min(kMaxCounterSize,
                                  stored_reads_size_
                                    + (k_ + 2) * current_reads_storage->size()));

    stored_reads_size_ = 0;
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

typedef tsl::hopscotch_map<KMer, uint32_t, KMerHash> Counter;

std::vector<bool> filter_reads(const std::vector<std::string> &reads,
                               Counter *counter,
                               size_t k,
                               size_t noise_kmer_frequency,
                               bool verbose) {
    size_t num_distinct_kmers = counter->size();
    size_t num_reads = reads.size();
    size_t num_frequent_reads = reads.size();
    size_t counter_sum = 0;
    size_t filtered_sum = 0;

    std::vector<bool> filter(reads.size(), true);

    for (auto it = counter->begin(); it != counter->end(); ) {
        counter_sum += it->second;
        if (it->second <= noise_kmer_frequency) {
            it = counter->erase(it);
        } else {
            filtered_sum += it->second;
            ++it;
        }
    }
    counter->rehash(counter->size());

    std::vector<KMer> read_kmers;

    for (size_t j = 0; j < reads.size(); ++j) {
        if (reads[j].size() <= k) {
            filter[j] = false;
            continue;
        }

        sequence_to_kmers(reads[j], k, &read_kmers, {});

        // iterate through all non-dummy k-mers
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            // all frequent k-mers are kept in counter here
            if (!counter->count(read_kmers[i])) {
                num_frequent_reads--;
                filter[j] = false;
                break;
            }
        }
        read_kmers.clear();
    }

    if (verbose) {
        std::cout << "\nFiltering out the k-mers collected...\n";
        std::cout << "Total reads:     " << num_reads << "\n";
        std::cout << "Filtered reads:  " << num_frequent_reads << "\n";
        std::cout << "Distinct k-mers:   " << num_distinct_kmers << std::endl;
        std::cout << "Total k-mers:      " << counter_sum << "\n";
        std::cout << "Frequent k-mers:   " << filtered_sum << std::endl;
    }
    return filter;
}

void extract_frequent_kmers(std::vector<std::string> *frequent_reads,
                            std::vector<std::string> *reads,
                            Counter *counter,
                            size_t k,
                            std::vector<KMer> *kmers,
                            size_t *end_sorted,
                            const std::vector<TAlphabet> &suffix,
                            size_t noise_kmer_frequency,
                            size_t num_threads,
                            bool verbose,
                            std::mutex *mutex) {
    // filter out reads with noisy k-mers
    auto filter = filter_reads(*reads, counter, k, noise_kmer_frequency, verbose);

    counter->clear();
    counter->rehash(0);

    for (size_t i = 0; i < reads->size(); ++i) {
        if (filter[i])
            frequent_reads->emplace_back(std::move(reads->at(i)));
    }
    reads->clear();

    // extract k-mers from the reads filtered
    extract_kmers(
        [frequent_reads](CallbackRead callback) {
            for (auto &&read : *frequent_reads) {
                callback(std::move(read));
            }
            frequent_reads->clear();
        },
        k, kmers, end_sorted, suffix, num_threads, verbose, mutex
    );
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
    counter.rehash(kMaxCounterSize);

    std::vector<std::string> reads;
    std::vector<std::string> frequent_reads;
    std::vector<KMer> read_kmers;

    generate_reads([&](const std::string &read) {
        read_kmers.clear();
        sequence_to_kmers(read, k, &read_kmers, {});
        // consider only non-dummy k-mers
        if (read_kmers.size() <= 2)
            return;

        if (counter.size() + read_kmers.size() > kMaxCounterSize) {
            extract_frequent_kmers(&frequent_reads, &reads, &counter, k, kmers, end_sorted, suffix,
                                   noise_kmer_frequency,
                                   num_threads, verbose, mutex);
            counter.rehash(kMaxCounterSize);
        }

        bool frequent = true;
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            if (++counter[read_kmers[i]] <= noise_kmer_frequency)
                frequent = false;
        }
        if (frequent) {
            frequent_reads.push_back(read);
        } else {
            reads.push_back(read);
        }
    });

    extract_frequent_kmers(&frequent_reads, &reads, &counter, k, kmers, end_sorted, suffix,
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


std::vector<bool> filter_reads(std::function<void(CallbackRead)> generate_reads,
                               size_t k,
                               size_t noise_kmer_frequency,
                               bool verbose) {
    std::vector<bool> filter;

    Counter counter;
    counter.rehash(kMaxCounterSize);

    std::vector<std::string> reads;
    std::vector<KMer> read_kmers;

    generate_reads([&](const std::string &read) {
        read_kmers.clear();
        sequence_to_kmers(read, k, &read_kmers, {});
        // consider only non-dummy k-mers
        if (read_kmers.size() <= 2)
            return;

        if (counter.size() + read_kmers.size() > kMaxCounterSize) {
            auto filter_block = filter_reads(reads, &counter, k,
                                             noise_kmer_frequency, verbose);
            reads.clear();
            counter.clear();
            counter.rehash(kMaxCounterSize);

            filter.insert(filter.end(), filter_block.begin(), filter_block.end());
        }

        reads.push_back(read);
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            counter[read_kmers[i]]++;
        }
    });

    if (reads.size()) {
        auto filter_block = filter_reads(
            reads, &counter, k,
            noise_kmer_frequency * counter.size() / kMaxCounterSize,
            verbose
        );
        filter.insert(filter.end(), filter_block.begin(), filter_block.end());
    }
    return filter;
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
