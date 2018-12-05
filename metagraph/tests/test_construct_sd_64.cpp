#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_sd.hpp"
#include "dbg_construct.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

typedef KMerPacked<uint64_t, KmerExtractor2Bit::kLogSigma> KMER;
const int kMaxK = sizeof(KMER) * 8 / KmerExtractor2Bit::kLogSigma - 1;


TEST(Construct_SD_64, ConstructionNEAppendingSimplePath) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGSD constructed(&constructor);
        EXPECT_EQ(1u, constructed.num_nodes());

        DBGSD appended(k);
        ASSERT_EQ(appended.capacity(), appended.num_nodes()) << k;
        appended.add_sequence(std::string(100, 'A'));
        ASSERT_EQ(appended.capacity(), appended.num_nodes()) << k;

        ASSERT_NE(constructed, appended);
    }
}

TEST(Construct_SD_64, ConstructionNEAppendingTwoPaths) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(50, 'C') });
        DBGSD constructed(&constructor);
        EXPECT_EQ(2u, constructed.num_nodes());

        DBGSD appended(k);
        appended.add_sequence(std::string(100, 'A'));
        ASSERT_EQ(appended.capacity(), appended.num_nodes()) << k;
        appended.add_sequence(std::string(50, 'C'));
        ASSERT_EQ(appended.capacity(), appended.num_nodes()) << k;

        ASSERT_NE(constructed, appended);
    }
}

// TODO
/*
TEST(Construct_SD_64, ConstructionLowerCase) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGSDConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C') });
        DBGSD first(&constructor_first);

        DBGSDConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(100, 'a'),
                                           std::string(50, 'c') });
        DBGSD second(&constructor_second);

#if _DNA_CASE_SENSITIVE_GRAPH
        EXPECT_FALSE(first.equals_internally(second));
#else
        EXPECT_TRUE(first.equals_internally(second));
#endif
    }
}
*/

TEST(Construct_SD_64, ConstructionEQAppending) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        DBGSDConstructor constructor(k);
        constructor.add_sequences(input_data);
        DBGSD constructed(&constructor);

        Vector<KmerExtractor2Bit::Kmer64> kmers;
        for (const auto &str : input_data) {
            KmerExtractor2Bit::sequence_to_kmers<KmerExtractor2Bit::Kmer64>(
                str, k, {}, &kmers
            );
        }
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        EXPECT_EQ(kmers.size(), constructed.num_nodes());

        DBGSD appended(k);
        for (const auto &sequence : input_data) {
            ASSERT_EQ(appended.capacity(), appended.num_nodes()) << k;
            appended.add_sequence(sequence);
        }

        EXPECT_NE(constructed, appended);
    }
}

using TAlphabet = KmerExtractor2Bit::TAlphabet;

TEST(ExtractKmersPacked_64, ExtractKmersEmptySuffix) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, {}, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 700; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, {}, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < 500; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ 3 });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, suffix, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersWithFilteringTwo) {
    std::vector<TAlphabet> suffix = { 3, 1 };
    for (size_t k = 3; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::vector<TAlphabet>(length, 3), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();

            std::vector<TAlphabet> sequence(length, 3);
            sequence[k - 1] = 1;

            KmerExtractor2Bit::sequence_to_kmers(
                sequence, k, suffix, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersAppend) {
    Vector<KMER> result;

    KmerExtractor2Bit::sequence_to_kmers(
        std::vector<TAlphabet>(500, 3), 2, {}, &result
    );
    ASSERT_EQ(499u, result.size());

    KmerExtractor2Bit::sequence_to_kmers(
        std::vector<TAlphabet>(500, 3), 2, {}, &result
    );
    ASSERT_EQ(499u * 2, result.size());
}

TEST(ExtractKmersPacked_64, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        // NNN -> $NNN$
        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, {}, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersFromStringWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ KmerExtractor2Bit::encode('T') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}


TEST(ExtractKmersPacked_64, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 3; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractor2Bit::encode("TT"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, { 0, KmerExtractor2Bit::encode('T') }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractor2Bit::encode("TT"), &result
            );

            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor2Bit::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractor2Bit::encode("TA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'T');
            sequence[k - 1] = 'A';

            KmerExtractor2Bit::sequence_to_kmers(
                sequence, k, KmerExtractor2Bit::encode("TA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersFromStringAppend) {
    Vector<KMER> result;

    KmerExtractor2Bit::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u, result.size());

    KmerExtractor2Bit::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u * 2, result.size());
}


typedef std::function<void(const std::string&)> CallbackString;

template <typename KMER, class KmerExtractor2Bit>
void extract_kmers(std::function<void(CallbackString)> generate_reads,
                   size_t k,
                   bool canonical_mode,
                   Vector<KMER> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex &mutex_resize,
                   std::shared_timed_mutex &mutex_copy,
                   bool remove_redundant = true);

// TODO: k is node length
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        Vector<KMER> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        std::mutex &mutex_resize,
                                        std::shared_timed_mutex &mutex,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->reserve(reserved_capacity);
    extract_kmers<KMER, KmerExtractor2Bit>(
        [reads](CallbackString callback) {
            for (auto &&read : *reads) {
                callback(std::move(read));
            }
        },
        k, false, kmers, suffix,
        1, false, std::ref(mutex_resize), std::ref(mutex), remove_redundant
    );
    delete reads;
}

template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads = 1,
                                size_t offset = 0);

TEST(ExtractKmersPacked_64, ExtractKmersAppendParallelReserved) {
    Vector<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size - 1) * 5, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size - 1) * 10, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size - 1) * 15, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size - 1) * 20, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size - 1) * 20, result.size());
}

TEST(ExtractKmersPacked_64, ExtractKmersAppendParallel) {
    Vector<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(1u, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    // B->A in 2Bit mode
    ASSERT_EQ(1u, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(2u, result.size());
}

TEST(ExtractKmersPacked_64, ExtractKmersParallelRemoveRedundantReserved) {
    Vector<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 100'000
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 100'000
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, mu_resize, mu, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, mu_resize, mu, true, 100'000
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, mu_resize, mu, true, 100'000
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, mu_resize, mu, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, mu_resize, mu, true, 100'000
    );
    ASSERT_EQ(1u, result.size());
}

TEST(ExtractKmersPacked_64, ExtractKmersParallelRemoveRedundant) {
    Vector<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 0
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 0
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, mu_resize, mu, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, mu_resize, mu, true, 0
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, mu_resize, mu, true, 0
    );
    ASSERT_EQ(1u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, mu_resize, mu, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, mu_resize, mu, true, 0
    );
    ASSERT_EQ(1u, result.size());
}
