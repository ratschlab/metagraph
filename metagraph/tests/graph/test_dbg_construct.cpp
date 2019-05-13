#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "boss.hpp"
#include "boss_merge.hpp"
#include "boss_construct.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

template <typename Kmer>
class Construct : public ::testing::Test { };

template <typename Kmer>
class ExtractKmers : public ::testing::Test { };

typedef ::testing::Types<KMerBOSS<uint64_t, KmerExtractor::kLogSigma>,
                         KMerBOSS<sdsl::uint128_t, KmerExtractor::kLogSigma>,
                         KMerBOSS<sdsl::uint256_t, KmerExtractor::kLogSigma>> KmerTypes;

TYPED_TEST_CASE(Construct, KmerTypes);
TYPED_TEST_CASE(ExtractKmers, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractor::kLogSigma )


TYPED_TEST(Construct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(Construct, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(50, 'B') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));
        appended.add_sequence(std::string(50, 'B'));

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(Construct, ConstructionLowerCase) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(100, 'a'),
                                           std::string(50, 'c') });
        BOSS second(&constructor_second);

#if _DNA_CASE_SENSITIVE_GRAPH
        EXPECT_FALSE(first.equals_internally(second));
#else
        EXPECT_TRUE(first.equals_internally(second));
#endif
    }
}

TYPED_TEST(Construct, ConstructionDummySentinel) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'N'),
                                          std::string(50, '$') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(100, 'N'),
                                           std::string(50, '.') });
        BOSS second(&constructor_second);

        EXPECT_TRUE(first.equals_internally(second));
    }
}

TYPED_TEST(Construct, ConstructionEQAppending) {
    for (size_t k = 1; k < kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        BOSSConstructor constructor(k);
        constructor.add_sequences(input_data);
        BOSS constructed(&constructor);

        BOSS appended(k);
        for (const auto &sequence : input_data) {
            appended.add_sequence(sequence);
        }

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(Construct, ConstructionEQAppendingCanonical) {
    for (size_t k = 1; k < kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        BOSSConstructor constructor(k, true);
        constructor.add_sequences(input_data);
        BOSS constructed(&constructor);

        BOSS appended(k);
        for (auto &sequence : input_data) {
            appended.add_sequence(sequence);
            reverse_complement(sequence.begin(), sequence.end());
            appended.add_sequence(sequence);
        }

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(Construct, ConstructionLong) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(k + 1, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(k + 1, 'A'));

        EXPECT_EQ(constructed, appended);
        ASSERT_TRUE(constructed.num_nodes() > 1u);
    }
}

TYPED_TEST(Construct, ConstructionShort) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(k, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(k, 'A'));

        EXPECT_EQ(constructed, appended);
        ASSERT_EQ(1u, constructed.num_nodes());
    }
}

TYPED_TEST(Construct, ConstructionFromChunks) {
    for (size_t k = 1; k < kMaxK; k += 6) {
        BOSS boss_dynamic(k);
        boss_dynamic.add_sequence(std::string(100, 'A'));
        boss_dynamic.add_sequence(std::string(100, 'C'));
        boss_dynamic.add_sequence(std::string(100, 'T') + "A"
                                        + std::string(100, 'G'));

        for (size_t suffix_len = 0; suffix_len < k && suffix_len <= 5u; ++suffix_len) {
            BOSS::Chunk graph_data(k);

            for (const std::string &suffix : KmerExtractor::generate_suffixes(suffix_len)) {
                std::unique_ptr<IBOSSChunkConstructor> constructor(
                    IBOSSChunkConstructor::initialize(k, false, suffix)
                );

                constructor->add_sequence(std::string(100, 'A'));
                constructor->add_sequence(std::string(100, 'C'));
                constructor->add_sequence(std::string(100, 'T') + "A"
                                                + std::string(100, 'G'));

                auto next_block = constructor->build_chunk();
                graph_data.extend(*next_block);
                delete next_block;
            }

            BOSS boss;
            graph_data.initialize_boss(&boss);

            EXPECT_EQ(boss_dynamic, boss);
        }
    }
}


using TAlphabet = KmerExtractor::TAlphabet;

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        // NNN -> $NNN$
        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, {}, &result
            );
            ASSERT_EQ(length - k + 3, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }

    suffix.assign({ KmerExtractor::encode('N') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        // NN   -> $NN$
        // NNN  -> $$NNN$
        // NNNN -> $$$NNNN$
        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_EQ(length, result.size()) << "k: " << k
                                             << ", length: " << length;
        }
    }
}


TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 3; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NN"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        // NNN -> $NNN$
        for (size_t length = k; length < 200; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, { 0, KmerExtractor::encode('N') }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NN"), &result
            );

            // NNN  -> $$NNN$
            // NNNN -> $$$NNNN$
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'N');
            sequence[k - 1] = 'A';

            KmerExtractor::sequence_to_kmers(
                sequence, k, KmerExtractor::encode("NA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringAppend) {
    Vector<TypeParam> result;

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


typedef std::function<void(const std::string&)> CallbackString;

template <typename TypeParam, class KmerExtractor>
void extract_kmers(std::function<void(CallbackString)> generate_reads,
                   size_t k,
                   bool canonical_mode,
                   Vector<TypeParam> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex &mutex_resize,
                   std::shared_timed_mutex &mutex_copy,
                   bool remove_redundant = true);

// TODO: k is node length
template <typename TypeParam>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        Vector<TypeParam> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        std::mutex &mutex_resize,
                                        std::shared_timed_mutex &mutex,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->reserve(reserved_capacity);
    extract_kmers<TypeParam, KmerExtractor>(
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

TYPED_TEST(ExtractKmers, ExtractKmersAppendParallelReserved) {
    Vector<TypeParam> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 5, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 10, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 15, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, mu_resize, mu, false, 100'000
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());
}

TYPED_TEST(ExtractKmers, ExtractKmersAppendParallel) {
    Vector<TypeParam> result;
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
    ASSERT_EQ(3u, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(6u, result.size());
}

TYPED_TEST(ExtractKmers, ExtractKmersParallelRemoveRedundantReserved) {
    Vector<TypeParam> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 100'000
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 100'000
    );
    // $AA, AAA, AA$
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
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, mu_resize, mu, true, 100'000
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, mu_resize, mu, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, mu_resize, mu, true, 100'000
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.size());
}

TYPED_TEST(ExtractKmers, ExtractKmersParallelRemoveRedundant) {
    Vector<TypeParam> result;
    std::mutex mu_resize;
    std::shared_timed_mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, mu_resize, mu, true, 0
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, mu_resize, mu, true, 0
    );
    // $AA, AAA, AA$
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
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, mu_resize, mu, true, 0
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, mu_resize, mu, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, mu_resize, mu, true, 0
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.size());
}
