#include "reads_filtering.hpp"

#include <numeric>

#include <tsl/hopscotch_map.h>
#include <kmc_file.h>

#include "kmer.hpp"
#include "utils.hpp"
#include "unix_tools.hpp"

using utils::KmerExtractor;

const size_t kMaxKmersChunkSize = 30'000'000;

template <typename KMER>
using Counter
    = typename tsl::hopscotch_pg_map<const KMER, uint32_t, utils::Hash<KMER>>;


template <typename KMER>
std::vector<bool> filter_reads(const std::vector<std::string> &reads,
                               Counter<KMER> *counter,
                               size_t k,
                               size_t max_unreliable_abundance,
                               size_t unreliable_kmers_threshold,
                               bool verbose) {
    Timer timer;

    size_t num_distinct_kmers = counter->size();
    size_t num_frequent_kmers = 0;
    size_t counter_sum = 0;
    size_t filtered_sum = 0;

    std::vector<bool> filter(reads.size(), true);

    for (auto it = counter->begin(); it != counter->end(); ) {
        counter_sum += it->second;
        if (it->second <= max_unreliable_abundance) {
            it = counter->erase(it);
        } else {
            filtered_sum += it->second;
            ++it;
            num_frequent_kmers++;
        }
    }
    counter->rehash(counter->size());

    Vector<KMER> read_kmers;

    for (size_t j = 0; j < reads.size(); ++j) {
        if (reads[j].size() <= k) {
            filter[j] = false;
            continue;
        }

        KmerExtractor::sequence_to_kmers(reads[j], k + 1, {}, &read_kmers);

        size_t unreliable_kmers = 0;

        // iterate through all non-dummy k-mers
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            // all frequent k-mers are kept in counter here
            if (!counter->count(read_kmers[i])
                    && ++unreliable_kmers > unreliable_kmers_threshold) {
                filter[j] = false;
                break;
            }
        }
        read_kmers.resize(0);
    }

    if (verbose) {
        std::cout << "Filtering out the k-mers collected... "
                                                << timer.elapsed() << "sec\n";
        std::cout << "Distinct k-mers:   " << num_distinct_kmers << "\n";
        std::cout << "Frequent k-mers:   " << num_frequent_kmers << "\n";
        std::cout << "Total k-mers:      " << counter_sum << "\n";
        std::cout << "Total k-mers left: " << filtered_sum << std::endl;
    }
    return filter;
}


// Identify k-mers from |reads| below the abundance level
// and return vector indicating reads with reliable k-mers
template <typename KMER>
std::vector<bool>
count_kmers_and_filter_reads_templated(std::vector<std::string> *reads,
                                       size_t k,
                                       size_t max_unreliable_abundance,
                                       size_t unreliable_kmers_threshold,
                                       bool verbose) {
    std::vector<bool> frequent(reads->size(), true);

    std::vector<std::string> filtering_reads;
    filtering_reads.reserve(reads->size());

    std::vector<std::string> frequent_reads;
    frequent_reads.reserve(reads->size());

    Counter<KMER> counter(std::accumulate(reads->begin(), reads->end(), 0,
        [&](size_t sum, const std::string &read) {
            return sum + std::max(read.size(), k) - k + 2;
        })
    );

    Vector<KMER> read_kmers;

    for (size_t j = 0; j < reads->size(); ++j) {
        KmerExtractor::sequence_to_kmers(reads->at(j), k + 1, {}, &read_kmers);

        // filter out too short reads
        if (read_kmers.size() <= 2)
            frequent[j] = false;

        // consider only non-dummy k-mers
        for (size_t i = 1; i + 1 < read_kmers.size(); ++i) {
            if (++counter[read_kmers[i]] <= max_unreliable_abundance)
                frequent[j] = false;
        }
        if (frequent[j]) {
            frequent_reads.emplace_back(std::move(reads->at(j)));
        } else {
            filtering_reads.emplace_back(std::move(reads->at(j)));
        }
        read_kmers.resize(0);
    }
    reads->clear();

    auto filter = filter_reads(filtering_reads, &counter, k,
                               max_unreliable_abundance,
                               unreliable_kmers_threshold, verbose);

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

template <typename... Args>
std::vector<bool> count_kmers_and_filter_reads(std::vector<std::string> *reads,
                                               size_t k,
                                               Args&... args) {
    if ((k + 1) * KmerExtractor::kLogSigma <= 64) {
        return count_kmers_and_filter_reads_templated<KmerExtractor::Kmer64>(
            reads, k, args...
        );
    } else if ((k + 1) * KmerExtractor::kLogSigma <= 128) {
        return count_kmers_and_filter_reads_templated<KmerExtractor::Kmer128>(
            reads, k, args...
        );
    } else {
        return count_kmers_and_filter_reads_templated<KmerExtractor::Kmer256>(
            reads, k, args...
        );
    }
}


// RAII, releases reads when finished
std::vector<bool>
count_kmers_and_filter_reads_raii(std::vector<std::string> *reads,
                                  size_t k,
                                  size_t max_unreliable_abundance,
                                  size_t unreliable_kmers_threshold,
                                  CKMCFile *kmc_database,
                                  bool verbose) {
    size_t all_reads = reads->size();
    std::vector<bool> result;

    if (kmc_database) {
        for (const std::string &read : *reads) {
            if (read.size() <= k) {
                result.push_back(0);
                continue;
            }

            std::vector<uint32> counters;
            size_t valid_kmers = 0;

            kmc_database->GetCountersForRead(read, counters);

            for (auto counter : counters) {
                if (counter)
                    ++valid_kmers;
            }

            result.push_back(valid_kmers + unreliable_kmers_threshold >= read.size() - k);
        }

    } else {
        result = count_kmers_and_filter_reads(
            reads, k, max_unreliable_abundance, unreliable_kmers_threshold, verbose
        );
    }

    delete reads;

    if (verbose) {
        std::cout << "All reads:   " << all_reads << "\n";
        std::cout << "Reads left:  "
                  << std::accumulate(result.begin(), result.end(), 0llu)
                  << "\n" << std::endl;
    }

    return result;
}


std::vector<bool> filter_reads(std::function<void(CallbackRead)> generate_reads,
                               size_t k,
                               size_t max_unreliable_abundance,
                               size_t unreliable_kmers_threshold,
                               bool verbose,
                               utils::ThreadPool *thread_pool,
                               const std::string &kmc_base) {
    std::unique_ptr<CKMCFile> kmc_database;

    if (kmc_base.size()) {
        kmc_database.reset(new CKMCFile());

        if (!kmc_database->OpenForRA(kmc_base)) {
            std::cerr << "Error: Can't open KMC database " << kmc_base << std::endl;
            exit(1);
        }
        if (kmc_database->KmerLength() != k + 1) {
            std::cerr << "Error: Incompatible KMC database " << kmc_base
                      << " with k=" << kmc_database->KmerLength()
                      << " instead of required k=" << k + 1 << std::endl;
            exit(1);
        }
        if (kmc_database->GetMinCount() > max_unreliable_abundance + 1) {
            std::cerr << "Error: Incompatible KMC database " << kmc_base
                      << " built with min_count=" << kmc_database->GetMinCount() << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "KMC database loaded" << std::endl;
        }

        kmc_database->SetMinCount(max_unreliable_abundance + 1);
    }

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
                    reads, k, max_unreliable_abundance, unreliable_kmers_threshold,
                    kmc_database.get(), verbose
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
            max_unreliable_abundance * valid_kmers_collected / kMaxKmersChunkSize,
            unreliable_kmers_threshold,
            kmc_database.get(), verbose
        )
    );

    std::vector<bool> reads_filter;
    for (auto &future_filter_block : future_filters) {
        auto filter_block = future_filter_block.get();
        reads_filter.insert(reads_filter.end(), filter_block.begin(), filter_block.end());
    }

    return reads_filter;
}
