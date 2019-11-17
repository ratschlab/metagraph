#include "graph/succinct/boss_construct.hpp"

#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "graph/succinct/boss.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_multiset.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

template <typename Kmer>
class BOSSConstruct : public ::testing::Test { };

template <typename Kmer>
class WeightedBOSSConstruct : public ::testing::Test { };

template <typename Kmer>
class CollectKmers : public ::testing::Test { };

template <typename Kmer>
class CountKmers : public ::testing::Test { };

template <typename KMER, bool Weighted = false>
class BOSSConfigurationType {
 public:
    typedef KMER Kmer;
    static const bool kWeighted = Weighted;
};
template <typename KMER, bool Weighted>
const bool BOSSConfigurationType<KMER, Weighted>::kWeighted;

typedef ::testing::Types<BOSSConfigurationType<KmerExtractorBOSS::Kmer64, false>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer128, false>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer256, false>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer64, true>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer128, true>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer256, true>> KmerAndWeightedTypes;

typedef ::testing::Types<BOSSConfigurationType<KmerExtractorBOSS::Kmer64, true>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer128, true>,
                         BOSSConfigurationType<KmerExtractorBOSS::Kmer256, true>> KmerWeightedTypes;

TYPED_TEST_CASE(BOSSConstruct, KmerAndWeightedTypes);
TYPED_TEST_CASE(WeightedBOSSConstruct, KmerWeightedTypes);

typedef ::testing::Types<KMerBOSS<uint64_t, KmerExtractorBOSS::bits_per_char>,
                         KMerBOSS<sdsl::uint128_t, KmerExtractorBOSS::bits_per_char>,
                         KMerBOSS<sdsl::uint256_t, KmerExtractorBOSS::bits_per_char>> KmerTypes;

TYPED_TEST_CASE(CollectKmers, KmerTypes);
TYPED_TEST_CASE(CountKmers, KmerTypes);

#define kMaxK ( sizeof(typename TypeParam::Kmer) * 8 / KmerExtractorBOSS::bits_per_char )


TYPED_TEST(BOSSConstruct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(BOSSConstruct, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
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
        BOSSConstructor constructor_first(k, false, TypeParam::kWeighted);
        constructor_first.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k, false, TypeParam::kWeighted);
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
        BOSSConstructor constructor_first(k, false, TypeParam::kWeighted);
        constructor_first.add_sequences({ std::string(100, 'N'),
                                          std::string(50, '$') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k, false, TypeParam::kWeighted);
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
        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
        constructor.add_sequences(input_data);
        BOSS constructed(&constructor);

        BOSS appended(k);
        for (const auto &sequence : input_data) {
            appended.add_sequence(sequence);
        }

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeight) {
    ASSERT_TRUE(TypeParam::kWeighted);

    for (size_t k = 1; k < kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };

        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
        constructor.add_sequences(input_data);

        BOSS constructed;
        sdsl::int_vector<> weights;
        constructor.build_graph(&constructed, &weights);

        ASSERT_EQ(constructed.num_edges() + 1, weights.size());

        auto mask = constructed.mark_all_dummy_edges(1);
        ASSERT_EQ(weights.size(), mask.size());

        for (size_t i = 1; i < weights.size(); ++i) {
            auto node_str = constructed.get_node_str(i)
                + constructed.decode(constructed.get_W(i) % constructed.alph_size);

            ASSERT_EQ(k + 1, node_str.size());
            ASSERT_EQ(node_str[0] == '$' || node_str[k] == '$', mask[i]);

            EXPECT_EQ(mask[i], weights[i] == 0) << i << " " << node_str;
        }
    }
}

TYPED_TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeightChunks) {
    ASSERT_TRUE(TypeParam::kWeighted);

    for (size_t k = 1; k < kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };

        BOSS constructed(k);

        auto constructor = IBOSSChunkConstructor::initialize(
            k,
            false,
            TypeParam::kWeighted
        );

        for (auto&& sequence : input_data) {
            constructor->add_sequence(std::move(sequence));
        }

        std::unique_ptr<BOSS::Chunk> chunk { constructor->build_chunk() };

        sdsl::int_vector<> weights;
        chunk->initialize_boss(&constructed, &weights);

        ASSERT_EQ(constructed.num_edges() + 1, weights.size());

        auto mask = constructed.mark_all_dummy_edges(1);
        ASSERT_EQ(weights.size(), mask.size());

        for (size_t i = 1; i < weights.size(); ++i) {
            auto node_str = constructed.get_node_str(i)
                + constructed.decode(constructed.get_W(i) % constructed.alph_size);

            ASSERT_EQ(k + 1, node_str.size());
            ASSERT_EQ(node_str[0] == '$' || node_str[k] == '$', mask[i]);

            EXPECT_EQ(mask[i], weights[i] == 0) << i << " " << node_str;
        }
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
        BOSSConstructor constructor(k, true, TypeParam::kWeighted);
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
        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
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
        BOSSConstructor constructor(k, false, TypeParam::kWeighted);
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
            BOSS::Chunk graph_data(KmerExtractorBOSS::alphabet.size(), k, false);

            for (const std::string &suffix : KmerExtractorBOSS::generate_suffixes(suffix_len)) {
                std::unique_ptr<IBOSSChunkConstructor> constructor(
                    IBOSSChunkConstructor::initialize(k, false, TypeParam::kWeighted, suffix)
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

TYPED_TEST(BOSSConstruct, ConstructionFromChunksParallel) {
    const uint64_t num_threads = 4;

    for (size_t k = 1; k < kMaxK; k += 6) {
        BOSS boss_dynamic(k);
        boss_dynamic.add_sequence(std::string(100, 'A'));
        boss_dynamic.add_sequence(std::string(100, 'C'));
        boss_dynamic.add_sequence(std::string(100, 'T') + "A"
                                        + std::string(100, 'G'));

        for (size_t suffix_len = 0; suffix_len < k && suffix_len <= 5u; ++suffix_len) {
            BOSS::Chunk graph_data(KmerExtractorBOSS::alphabet.size(), k, false);

            for (const std::string &suffix : KmerExtractorBOSS::generate_suffixes(suffix_len)) {
                std::unique_ptr<IBOSSChunkConstructor> constructor(
                    IBOSSChunkConstructor::initialize(k, false, TypeParam::kWeighted, suffix, num_threads)
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


template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix,
                   bool remove_redundant = true);

// TODO: k is node length
template <typename KMER>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        SortedSet<KMER, Vector<KMER>> *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    extract_kmers<KMER, KmerExtractorBOSS, SortedSet<KMER, Vector<KMER>>>(
        [reads](CallString callback) {
            std::for_each(reads->begin(), reads->end(), callback);
        },
        k, false, kmers, suffix, remove_redundant
    );
    delete reads;
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallelReserved) {
    SortedSet<TypeParam, Vector<TypeParam>> result;
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
#if _DNA_GRAPH
    ASSERT_EQ(3u, result.data().size());
#else
    ASSERT_EQ(6u, result.data().size());
#endif

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, false, 100'000
    );
#if _DNA_GRAPH
    ASSERT_EQ(3u, result.data().size());
#else
    ASSERT_EQ(6u, result.data().size());
#endif
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallel) {
    SortedSet<TypeParam, Vector<TypeParam>> result;
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
#if _DNA_GRAPH
    ASSERT_EQ(3u, result.data().size());
#else
    ASSERT_EQ(6u, result.data().size());
#endif
}

TYPED_TEST(CollectKmers, CollectKmersParallelRemoveRedundantReserved) {
    SortedSet<TypeParam, Vector<TypeParam>> result;

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
    SortedSet<TypeParam, Vector<TypeParam>> result;

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

typedef std::function<void(const std::string&, uint64_t)> CallStringCount;

template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix);

// TODO: k is node length
template <typename TypeParam, typename KmerCount>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        SortedMultiset<TypeParam,
                                                       KmerCount,
                                                       Vector<std::pair<TypeParam, KmerCount>>> *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    count_kmers<TypeParam,
                KmerExtractorBOSS,
                SortedMultiset<TypeParam,
                               KmerCount,
                               Vector<std::pair<TypeParam, KmerCount>>>>(
        [reads](CallStringCount callback) {
            for (const auto &read : *reads) {
                callback(read, 1);
            }
        },
        k, false, kmers, suffix
    );
    delete reads;
}

TYPED_TEST(CountKmers, CountKmers8bits) {
    SortedMultiset<TypeParam, uint8_t, Vector<std::pair<TypeParam, uint8_t>>> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(5u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(255u, result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(10u, result.data()[0].second);
    EXPECT_EQ(10u, result.data()[1].second);
    EXPECT_EQ(255u, result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(15u, result.data()[1].second);
    EXPECT_EQ(255u, result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(255u, result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(255u, result.data()[5].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(255u, result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(255u, result.data()[5].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 0, }, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(10u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(255u, result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(255u, result.data()[5].second);
}

TYPED_TEST(CountKmers, CountKmers32bits) {
    SortedMultiset<TypeParam, uint32_t, Vector<std::pair<TypeParam, uint32_t>>> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(5u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(5 * (sequence_size - 2 + 1), result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(10u, result.data()[0].second);
    EXPECT_EQ(10u, result.data()[1].second);
    EXPECT_EQ(10 * (sequence_size - 2 + 1), result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(3u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(15u, result.data()[1].second);
    EXPECT_EQ(15 * (sequence_size - 2 + 1), result.data()[2].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(15 * (sequence_size - 2 + 1), result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(5 * (sequence_size - 2 + 1), result.data()[5].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(5u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(15 * (sequence_size - 2 + 1), result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(5 * (sequence_size - 2 + 1), result.data()[5].second);

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 0, }, 100'000
    );
    ASSERT_EQ(6u, result.data().size());

    EXPECT_EQ(15u, result.data()[0].second);
    EXPECT_EQ(10u, result.data()[1].second);
    EXPECT_EQ(15u, result.data()[2].second);
    EXPECT_EQ(15 * (sequence_size - 2 + 1), result.data()[3].second);
    EXPECT_EQ(5u, result.data()[4].second);
    EXPECT_EQ(5 * (sequence_size - 2 + 1), result.data()[5].second);
}

TYPED_TEST(CountKmers, CountKmersAppendParallel) {
    SortedMultiset<TypeParam, uint8_t, Vector<std::pair<TypeParam, uint8_t>>> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, 0
    );
#if _DNA_GRAPH
    ASSERT_EQ(3u, result.data().size());
#else
    ASSERT_EQ(6u, result.data().size());
#endif
}
