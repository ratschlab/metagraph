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

typedef KMer<sdsl::uint256_t, KmerExtractor::kLogSigma> KMER;
const int kMaxK = sizeof(KMER) * 8 / KmerExtractor::kLogSigma;


TEST(Construct_256, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TEST(Construct_256, ConstructionEQAppendingTwoPaths) {
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

TEST(Construct_256, ConstructionLowerCase) {
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

TEST(Construct_256, ConstructionDummySentinel) {
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

TEST(Construct_256, ConstructionEQAppending) {
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

TEST(Construct_256, ConstructionLong) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(k + 1, 'A') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(k + 1, 'A'));

        EXPECT_EQ(constructed, appended);
        ASSERT_TRUE(constructed.num_nodes() > 1u);
    }
}

TEST(Construct_256, ConstructionShort) {
    for (size_t k = 1; k < kMaxK; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(k, 'A') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(k, 'A'));

        EXPECT_EQ(constructed, appended);
        ASSERT_EQ(1u, constructed.num_nodes());
    }
}

using TAlphabet = KmerExtractor::TAlphabet;

TEST(ExtractKmers_256, ExtractKmersEmptySuffix) {
    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, {}, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 700; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, {}, &result
            );
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_256, ExtractKmersWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < 500; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ 6 });
    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, suffix, &result
            );
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_256, ExtractKmersWithFilteringTwo) {
    std::vector<TAlphabet> suffix = { 6, 1 };
    for (size_t k = 2; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length <= k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::vector<TAlphabet>(length, 6), k + 1, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 500; ++length) {
            result.clear();

            std::vector<TAlphabet> sequence(length, 6);
            sequence[k - 1] = 1;

            KmerExtractor::sequence_to_kmers(
                sequence, k + 1, suffix, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_256, ExtractKmersAppend) {
    Vector<KMER> result;

    KmerExtractor::sequence_to_kmers(
        std::vector<TAlphabet>(500, 6), 2, {}, &result
    );
    ASSERT_EQ(499u, result.size());

    KmerExtractor::sequence_to_kmers(
        std::vector<TAlphabet>(500, 6), 2, {}, &result
    );
    ASSERT_EQ(499u * 2, result.size());
}

TEST(ExtractKmers_256, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 1; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length <= k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k + 1; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, {}, &result
            );
            // NNN -> $NNN$
            ASSERT_EQ(length - k + 2, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_256, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 2; k < kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length <= k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, KmerExtractor::encode("NN"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k + 1; length < 200; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, { 0, 0 }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, { 0, KmerExtractor::encode('N') }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, KmerExtractor::encode("NN"), &result
            );
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length <= k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k + 1, KmerExtractor::encode("NA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'N');
            sequence[k - 1] = 'A';

            KmerExtractor::sequence_to_kmers(
                sequence, k + 1, KmerExtractor::encode("NA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers_256, ExtractKmersFromStringAppend) {
    Vector<KMER> result;

    // A...A -> $A...A$
    KmerExtractor::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(501u, result.size());

    KmerExtractor::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(501u * 2, result.size());
}


typedef std::function<void(const std::string&)> CallbackRead;

template <typename KMER>
void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   Vector<KMER> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   const KmerExtractor &kmer_extractor,
                   size_t num_threads,
                   bool verbose,
                   std::mutex &mutex_resize,
                   std::shared_timed_mutex &mutex_copy,
                   bool remove_redundant = true);

void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        Vector<KMER> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        std::mutex &mutex_resize,
                                        std::shared_timed_mutex &mutex,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->reserve(reserved_capacity);
    KmerExtractor kmer_extractor;
    extract_kmers([reads](CallbackRead callback) {
            for (auto &&read : *reads) {
                callback(std::move(read));
            }
        },
        k, kmers, suffix, kmer_extractor,
        1, false, std::ref(mutex_resize), std::ref(mutex), remove_redundant
    );
    delete reads;
}

template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads = 1,
                                size_t offset = 0);

TEST(ExtractKmers_256, ExtractKmersAppendParallelReserved) {
    Vector<KMER> result;
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

TEST(ExtractKmers_256, ExtractKmersAppendParallel) {
    Vector<KMER> result;
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

TEST(ExtractKmers_256, ExtractKmersParallelRemoveRedundantReserved) {
    Vector<KMER> result;
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

TEST(ExtractKmers_256, ExtractKmersParallelRemoveRedundant) {
    Vector<KMER> result;
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
