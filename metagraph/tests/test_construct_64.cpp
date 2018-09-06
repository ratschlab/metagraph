#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "dbg_succinct_merge.hpp"
#include "dbg_succinct_construct.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

typedef KMer<uint64_t> KMER;
const int kMaxK = sizeof(KMER) * 8 / kBitsPerChar;


TEST(Construct_64, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TEST(Construct_64, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A'),
                                std::string(50, 'B') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(100, 'A'));
        appended.add_sequence(std::string(50, 'B'));

        EXPECT_EQ(constructed, appended);
    }
}

TEST(Construct_64, ConstructionLowerCase) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor_first(k);
        constructor_first.add_reads({ std::string(100, 'A'),
                                      std::string(50, 'C') });
        DBG_succ first(&constructor_first);

        KMerDBGSuccConstructor constructor_second(k);
        constructor_second.add_reads({ std::string(100, 'a'),
                                       std::string(50, 'c') });
        DBG_succ second(&constructor_second);

#if _DNA_CASE_SENSITIVE_GRAPH
        EXPECT_FALSE(first.equals_internally(second));
#else
        EXPECT_TRUE(first.equals_internally(second));
#endif
    }
}

TEST(Construct_64, ConstructionDummySentinel) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor_first(k);
        constructor_first.add_reads({ std::string(100, 'N'),
                                      std::string(50, '$') });
        DBG_succ first(&constructor_first);

        KMerDBGSuccConstructor constructor_second(k);
        constructor_second.add_reads({ std::string(100, 'N'),
                                       std::string(50, '.') });
        DBG_succ second(&constructor_second);

        EXPECT_TRUE(first.equals_internally(second));
    }
}

TEST(Construct_64, ConstructionEQAppending) {
    for (size_t k = 1; k < kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads(input_data);
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        for (const auto &sequence : input_data) {
            appended.add_sequence(sequence);
        }

        EXPECT_EQ(constructed, appended);
    }
}

namespace utils {
    template <typename KMER>
    void sequence_to_kmers(std::vector<TAlphabet>&& seq,
                           size_t k,
                           Vector<KMER> *kmers,
                           const std::vector<TAlphabet> &suffix);
}

TEST(ExtractKmers_64, ExtractKmersEmptySuffix) {
    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, {});
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 700; ++length) {
            result.clear();
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, {});
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_64, ExtractKmersWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < 500; ++length) {
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ 6 });
    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_64, ExtractKmersWithFilteringTwo) {
    std::vector<TAlphabet> suffix = { 6, 1 };
    for (size_t k = 2; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length <= k; ++length) {
            utils::sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 500; ++length) {
            result.clear();

            std::vector<TAlphabet> sequence(length, 6);
            sequence[k - 1] = 1;

            utils::sequence_to_kmers(std::move(sequence), k, &result, suffix);
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_64, ExtractKmersAppend) {
    Vector<KMER> result;

    utils::sequence_to_kmers(std::vector<TAlphabet>(500, 6), 1, &result, {});
    ASSERT_EQ(499u, result.size());

    utils::sequence_to_kmers(std::vector<TAlphabet>(500, 6), 1, &result, {});
    ASSERT_EQ(499u * 2, result.size());
}

TEST(ExtractKmers_64, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result, {});
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result, {});
            // NNN -> $NNN$
            ASSERT_EQ(length - k + 2, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_64, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 2; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result,
                                     { DBG_succ::encode('N'),
                                       DBG_succ::encode('N') });
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result, { 0, 0 });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result,
                                     { 0, DBG_succ::encode('N') });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result,
                                     { DBG_succ::encode('N'),
                                       DBG_succ::encode('N') });
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length <= k; ++length) {
            utils::sequence_to_kmers(std::string(length, 'N'), k, &result,
                                     { DBG_succ::encode('N'), 1 });
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'N');
            sequence[k - 1] = 'A';

            utils::sequence_to_kmers(sequence, k, &result,
                                     { DBG_succ::encode('N'), 1 });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_64, ExtractKmersFromStringAppend) {
    Vector<KMER> result;

    utils::sequence_to_kmers(std::string(500, 'A'), 1, &result, {});
    ASSERT_EQ(501u, result.size());

    utils::sequence_to_kmers(std::string(500, 'A'), 1, &result, {});
    ASSERT_EQ(501u * 2, result.size());
}


typedef std::function<void(const std::string&)> CallbackRead;

template <typename KMER>
void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   utils::DequeStorage<KMER> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex &mutex_resize,
                   std::shared_timed_mutex &mutex,
                   bool remove_redundant = true);

void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        utils::DequeStorage<KMER> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        std::mutex &mutex_resize,
                                        std::shared_timed_mutex &mutex,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->reserve(reserved_capacity);
    extract_kmers([reads](CallbackRead callback) {
            for (auto &&read : *reads) {
                callback(std::move(read));
            }
        },
        k, kmers, suffix,
        1, false, std::ref(mutex_resize), std::ref(mutex), remove_redundant
    );
    delete reads;
}

template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads = 1,
                                size_t offset = 0);

TEST(ExtractKmers_64, ExtractKmersAppendParallelReserved) {
    utils::DequeStorage<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 5, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 10, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 15, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, { 1, }, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());
}

TEST(ExtractKmers_64, ExtractKmersAppendParallel) {
    utils::DequeStorage<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(3u, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(6u, result.size());
}

TEST(ExtractKmers_64, ExtractKmersParallelRemoveRedundantReserved) {
    utils::DequeStorage<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        1, &result, {}, mu_resize, mu, true, 100'000
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 100'000
    );
    // $AA, AAA, AA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 0 }, mu_resize, mu, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 1 }, mu_resize, mu, true, 100'000
    );
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 100'000
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, mu_resize, mu, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, mu_resize, mu, true, 100'000
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.size());
}

TEST(ExtractKmers_64, ExtractKmersParallelRemoveRedundant) {
    utils::DequeStorage<KMER> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        1, &result, {}, mu_resize, mu, true, 0
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 0
    );
    // $AA, AAA, AA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 0 }, mu_resize, mu, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 1 }, mu_resize, mu, true, 0
    );
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 0
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, mu_resize, mu, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, mu_resize, mu, true, 0
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.size());
}
