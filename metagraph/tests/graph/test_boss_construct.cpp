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
#include "sorted_set.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

template <typename Kmer>
class BOSSConstruct : public ::testing::Test { };

template <typename Kmer>
class CollectKmers : public ::testing::Test { };

typedef ::testing::Types<KMerBOSS<uint64_t, KmerExtractor::kLogSigma>,
                         KMerBOSS<sdsl::uint128_t, KmerExtractor::kLogSigma>,
                         KMerBOSS<sdsl::uint256_t, KmerExtractor::kLogSigma>> KmerTypes;

TYPED_TEST_CASE(BOSSConstruct, KmerTypes);
TYPED_TEST_CASE(CollectKmers, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractor::kLogSigma )


TYPED_TEST(BOSSConstruct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(BOSSConstruct, ConstructionEQAppendingTwoPaths) {
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

TYPED_TEST(BOSSConstruct, ConstructionLowerCase) {
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

TYPED_TEST(BOSSConstruct, ConstructionDummySentinel) {
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

TYPED_TEST(BOSSConstruct, ConstructionEQAppending) {
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

TYPED_TEST(BOSSConstruct, ConstructionEQAppendingCanonical) {
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

TYPED_TEST(BOSSConstruct, ConstructionLong) {
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

TYPED_TEST(BOSSConstruct, ConstructionShort) {
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

TYPED_TEST(BOSSConstruct, ConstructionFromChunks) {
    for (size_t k = 1; k < kMaxK; k += 6) {
        BOSS boss_dynamic(k);
        boss_dynamic.add_sequence(std::string(100, 'A'));
        boss_dynamic.add_sequence(std::string(100, 'C'));
        boss_dynamic.add_sequence(std::string(100, 'T') + "A"
                                        + std::string(100, 'G'));

        for (size_t suffix_len = 0; suffix_len < k && suffix_len <= 5u; ++suffix_len) {
            BOSS::Chunk graph_data(KmerExtractor::alphabet.size(), k);

            for (const std::string &suffix : KmerExtractor::generate_suffixes(suffix_len)) {
                std::unique_ptr<IBOSSChunkConstructor> constructor(
                    IBOSSChunkConstructor::initialize(k, false, false, suffix)
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

typedef std::function<void(const std::string&)> CallbackString;

template <typename TypeParam, class KmerExtractor>
void extract_kmers(std::function<void(CallbackString)> generate_reads,
                   size_t k,
                   bool canonical_mode,
                   SortedSet<TypeParam> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   bool remove_redundant = true);

// TODO: k is node length
template <typename TypeParam>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        SortedSet<TypeParam> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    extract_kmers<TypeParam, KmerExtractor>(
        [reads](CallbackString callback) {
            std::for_each(reads->begin(), reads->end(), callback);
        },
        k, false, kmers, suffix, remove_redundant
    );
    delete reads;
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallelReserved) {
    SortedSet<TypeParam> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, false, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, false, 100'000
    );
    ASSERT_EQ(6u, result.data().size());
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallel) {
    SortedSet<TypeParam> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, false, 0
    );
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, false, 0
    );
    ASSERT_EQ(6u, result.data().size());
}

TYPED_TEST(CollectKmers, CollectKmersParallelRemoveRedundantReserved) {
    SortedSet<TypeParam> result;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, true, 100'000
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, true, 100'000
    );
    // $AA, AAA, AA$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, true, 100'000
    );
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, true, 100'000
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, true, 100'000
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, true, 100'000
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.data().size());
}

TYPED_TEST(CollectKmers, CollectKmersParallelRemoveRedundant) {
    SortedSet<TypeParam> result;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, true, 0
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, true, 0
    );
    // $AA, AAA, AA$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, true, 0
    );
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, {}, true, 0
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.data().size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 0 }, true, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        4, &result, { 1 }, true, 0
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.data().size());
}
