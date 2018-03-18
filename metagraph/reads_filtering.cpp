#include "reads_filtering.hpp"

#include <hopscotch_map.h>

#include "kmer.hpp"
#include "utils.hpp"


const size_t kMaxKmersChunkSize = 30'000'000;


typedef tsl::hopscotch_map<KMer, uint32_t, utils::KMerHash> Counter;

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

        utils::sequence_to_kmers(reads[j], k, &read_kmers, {});

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

// Remove noisy k-mers from |reads| and return
// vector indicating reads with frequent k-mers
std::vector<bool> count_kmers_and_filter_reads(std::vector<std::string> *reads,
                                               size_t k,
                                               size_t noise_kmer_frequency,
                                               bool verbose) {
    std::vector<bool> frequent(reads->size(), true);

    std::vector<std::string> filtering_reads;
    filtering_reads.reserve(reads->size());

    std::vector<std::string> frequent_reads;
    frequent_reads.reserve(reads->size());

    Counter counter;
    counter.rehash(kMaxKmersChunkSize / 3);

    std::vector<KMer> read_kmers;

    for (size_t j = 0; j < reads->size(); ++j) {
        utils::sequence_to_kmers(reads->at(j), k, &read_kmers, {});

        // consider only non-dummy k-mers
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            if (++counter[read_kmers[i]] <= noise_kmer_frequency)
                frequent[j] = false;
        }
        if (frequent[j]) {
            frequent_reads.emplace_back(std::move(reads->at(j)));
        } else {
            filtering_reads.emplace_back(std::move(reads->at(j)));
        }
        read_kmers.clear();
    }
    reads->clear();

    auto filter = filter_reads(filtering_reads, &counter,
                               k, noise_kmer_frequency, verbose);

    for (size_t i = 0, j = 0; j < frequent.size(); ++j) {
        if (!frequent[j] && filter[i++]) {
            frequent[j] = true;
            frequent_reads.emplace_back(std::move(filtering_reads[i - 1]));
        }
    }
    filtering_reads.clear();

    reads->swap(frequent_reads);
    return frequent;
}


// RAII, releases reads when finished
std::vector<bool>
count_kmers_and_filter_reads_raii(std::vector<std::string> *reads,
                                  size_t k,
                                  size_t noise_kmer_frequency,
                                  bool verbose) {
    auto result = count_kmers_and_filter_reads(
        reads, k, noise_kmer_frequency, verbose
    );
    delete reads;
    return result;
}

std::vector<bool> filter_reads(std::function<void(CallbackRead)> generate_reads,
                               size_t k,
                               size_t noise_kmer_frequency,
                               bool verbose,
                               utils::ThreadPool *thread_pool) {
    std::vector<std::future<std::vector<bool>>> future_filters;

    auto *reads = new std::vector<std::string>();
    size_t valid_kmers_collected = 0;

    generate_reads([&](const std::string &read) {
        size_t num_valid_kmers = read.size() > k
                                    ? read.size() - k
                                    : 0;

        if (valid_kmers_collected + num_valid_kmers > kMaxKmersChunkSize) {
            future_filters.emplace_back(
                thread_pool->enqueue(
                    count_kmers_and_filter_reads_raii,
                    reads, k, noise_kmer_frequency, verbose
                )
            );
            reads = new std::vector<std::string>();
            valid_kmers_collected = 0;
        }

        reads->push_back(read);
        valid_kmers_collected += num_valid_kmers;
    });

    future_filters.emplace_back(
        thread_pool->enqueue(
            count_kmers_and_filter_reads_raii,
            reads, k,
            noise_kmer_frequency * valid_kmers_collected / kMaxKmersChunkSize,
            verbose
        )
    );

    std::vector<bool> reads_filter;
    for (auto &future_filter_block : future_filters) {
        auto filter_block = future_filter_block.get();
        reads_filter.insert(reads_filter.end(), filter_block.begin(), filter_block.end());
    }

    return reads_filter;
}
