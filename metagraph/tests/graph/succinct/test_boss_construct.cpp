#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_multiset_disk.hpp"
#include "kmer/kmer_collector.hpp"

namespace {
using namespace mg;
using namespace mg::succinct;

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = TEST_DATA_DIR;
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

TYPED_TEST_SUITE(BOSSConstruct, KmerAndWeightedTypes);
TYPED_TEST_SUITE(WeightedBOSSConstruct, KmerWeightedTypes);

typedef ::testing::Types<KMerBOSS<uint64_t, KmerExtractorBOSS::bits_per_char>,
                         KMerBOSS<sdsl::uint128_t, KmerExtractorBOSS::bits_per_char>,
                         KMerBOSS<sdsl::uint256_t, KmerExtractorBOSS::bits_per_char>> KmerTypes;

TYPED_TEST_SUITE(CollectKmers, KmerTypes);
TYPED_TEST_SUITE(CountKmers, KmerTypes);

#define kMaxK ( sizeof(typename TypeParam::Kmer) * 8 / KmerExtractorBOSS::bits_per_char )


TYPED_TEST(BOSSConstruct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS constructed(&constructor);

        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TYPED_TEST(BOSSConstruct, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0);
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
        BOSSConstructor constructor_first(k, false, TypeParam::kWeighted ? 8 : 0);
        constructor_first.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k, false, TypeParam::kWeighted ? 8 : 0);
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
        BOSSConstructor constructor_first(k, false, TypeParam::kWeighted ? 8 : 0);
        constructor_first.add_sequences({ std::string(100, 'N'),
                                          std::string(50, '$') });
        BOSS first(&constructor_first);

        BOSSConstructor constructor_second(k, false, TypeParam::kWeighted ? 8 : 0);
        constructor_second.add_sequences({ std::string(100, 'N'),
                                           std::string(50, '.') });
        BOSS second(&constructor_second);

        EXPECT_TRUE(first.equals_internally(second));
    }
}

TYPED_TEST(BOSSConstruct, ConstructionEQAppending) {
    common::logger->set_level(spdlog::level::trace);
    for (auto container : { kmer::ContainerType::VECTOR_DISK }) { //TODO: undo change
        for (size_t k = 2; k < 3; ++k) {
            std::vector<std::string> input_data = {
                "ACAGCTAGCTAGCTAGCTAGCTG",
                "ATATTATAAAAAATTTTAAAAAA",
                "ATATATTCTCTCTCTCTCATA",
                "GTGTGTGTGGGGGGCCCTTTTTTCATA",
            };
            BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0, "", 1,
                                        20000, container);
            constructor.add_sequences(input_data);
            BOSS constructed(&constructor);

            BOSS appended(k);
            for (const auto &sequence : input_data) {
                appended.add_sequence(sequence);
            }

            EXPECT_EQ(constructed, appended);
        }
    }
}

TYPED_TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeight) {
    ASSERT_TRUE(TypeParam::kWeighted);
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            std::vector<std::string> input_data = {
                "ACAGCTAGCTAGCTAGCTAGCTG",
                "ATATTATAAAAAATTTTAAAAAA",
                "ATATATTCTCTCTCTCTCATA",
                "GTGTGTGTGGGGGGCCCTTTTTTCATA",
            };

            BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0, "", 1,
                                        20000, container);
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
}

TYPED_TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeightChunks) {
    ASSERT_TRUE(TypeParam::kWeighted);

    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            std::vector<std::string> input_data = {
                "ACAGCTAGCTAGCTAGCTAGCTG",
                "ATATTATAAAAAATTTTAAAAAA",
                "ATATATTCTCTCTCTCTCATA",
                "GTGTGTGTGGGGGGCCCTTTTTTCATA",
            };

            BOSS constructed(k);

            auto constructor
                    = IBOSSChunkConstructor::initialize(k, false, TypeParam::kWeighted ? 8 : 0,
                                                        "", 1, 20000, container);

            for (auto &&sequence : input_data) {
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
}

TYPED_TEST(BOSSConstruct, ConstructionEQAppendingCanonical) {
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            std::vector<std::string> input_data = {
                "ACAGCTAGCTAGCTAGCTAGCTG",
                "ATATTATAAAAAATTTTAAAAAA",
                "ATATATTCTCTCTCTCTCATA",
                "GTGTGTGTGGGGGGCCCTTTTTTCATA",
            };
            BOSSConstructor constructor(k, true, TypeParam::kWeighted ? 8 : 0, "", 1,
                                        20'000, container);
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
}

TYPED_TEST(BOSSConstruct, ConstructionLong) {
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0, "", 1,
                                        20'000, container);
            constructor.add_sequences({ std::string(k + 1, 'A') });
            BOSS constructed(&constructor);

            BOSS appended(k);
            appended.add_sequence(std::string(k + 1, 'A'));

            EXPECT_EQ(constructed, appended);
            ASSERT_TRUE(constructed.num_nodes() > 1u);
        }
    }
}

TYPED_TEST(BOSSConstruct, ConstructionShort) {
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            BOSSConstructor constructor(k, false, TypeParam::kWeighted ? 8 : 0, "", 1,
                                        20'000, container);
            constructor.add_sequences({ std::string(k, 'A') });
            BOSS constructed(&constructor);

            BOSS appended(k);
            appended.add_sequence(std::string(k, 'A'));

            EXPECT_EQ(constructed, appended);
        ASSERT_EQ(1u, constructed.num_nodes());
        }
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
                        IBOSSChunkConstructor::initialize(k, false, TypeParam::kWeighted ? 8 : 0,
                                                          suffix));

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
                        IBOSSChunkConstructor::initialize(k, false, TypeParam::kWeighted ? 8 : 0,
                                                          suffix, num_threads));

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


// TODO: k is node length
template <typename KMER>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        common::SortedSet<KMER, Vector<KMER>> *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix,
                                        bool remove_redundant,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    kmer::extract_kmers<KMER, KmerExtractorBOSS, common::SortedSet<KMER, Vector<KMER>>>(
        [reads](kmer::CallString callback) {
            std::for_each(reads->begin(), reads->end(), callback);
        },
        k, false, kmers, suffix, remove_redundant
    );
    delete reads;
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallelReserved) {
    common::SortedSet<TypeParam, Vector<TypeParam>> result;
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
    common::SortedSet<TypeParam, Vector<TypeParam>> result;
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
    common::SortedSet<TypeParam, Vector<TypeParam>> result;

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
    common::SortedSet<TypeParam, Vector<TypeParam>> result;

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

// TODO: k is node length
template <typename Container>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        Container *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix) {
    // kmers->try_reserve(reserved_capacity);
    kmer::count_kmers<typename Container::key_type, KmerExtractorBOSS, Container>(
            [reads](kmer::CallStringCount callback) {
            for (const auto &read : *reads) {
                callback(read, 1);
            }
        },
        k, false, kmers, suffix
    );
    delete reads;
}

/**
 * Checks that the given container contains the expected list of elements
 * @tparam Container the container to check the contents for, in our case a
 * SortedMultiset or a SortedMultisetDisk
 */
template <typename Container>
void assert_contents(Container &c, const std::initializer_list<size_t> &expected_els) {
    size_t size = 0;
    auto it_els = expected_els.begin();
    for (typename Container::result_type::iterator it = c.data().begin();
         it != c.data().end(); ++it, ++size, ++it_els) {
        EXPECT_EQ(*it_els, (*it).second);
    }
    EXPECT_EQ(expected_els.size(), size);
}

template <typename Container>
void check_counts() {
    std::function<void(typename Container::storage_type *)> cleanup
            = [](typename Container::storage_type *) {};
    Container result(cleanup, 1, 100'000);
    size_t sequence_size = 500;
    size_t max_value = std::numeric_limits<typename Container::count_type>::max();
    const size_t five_times = std::min(5 * (sequence_size - 2 + 1), max_value);
    const size_t ten_times = std::min(10 * (sequence_size - 2 + 1), max_value);
    const size_t fifteen_times = std::min(15 * (sequence_size - 2 + 1), max_value);

    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 5u, 5u, five_times });


    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 10u, 10u, ten_times });


    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 15u, 15u, fifteen_times });


    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result, {});
    assert_contents(result, { 15u, 5u, 15u, fifteen_times, 5u, five_times });


    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result,
            { 1 });
    assert_contents(result, { 15u, 5u, 15u, fifteen_times, 5u, five_times });


    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result,
            { 0 });
    assert_contents(result, { 15u, 10u, 15u, fifteen_times, 5u, five_times });
}

TYPED_TEST(CountKmers, CountKmers8bits) {
    using Container
            = common::SortedMultiset<TypeParam, uint8_t, Vector<std::pair<TypeParam, uint8_t>>>;
    check_counts<Container>();
}

TYPED_TEST(CountKmers, CountKmers32bits) {
    using Container
            = common::SortedMultiset<TypeParam, uint32_t, Vector<std::pair<TypeParam, uint32_t>>>;
    check_counts<Container>();
}

TYPED_TEST(CountKmers, CountKmers8bitsDisk) {
    using Container = common::SortedMultisetDisk<TypeParam, typename TypeParam::WordType, uint8_t>;
    std::function<void(typename Container::storage_type *)> cleanup
            = [](typename Container::storage_type *) {};
    Container result(cleanup, 1, 100'000);
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 5u, 5u, 255u });
}

TYPED_TEST(CountKmers, CountKmers32bitsDisk) {
    using Container = common::SortedMultisetDisk<TypeParam, typename TypeParam::WordType, uint32_t>;
    std::function<void(typename Container::storage_type *)> cleanup
            = [](typename Container::storage_type *) {};
    Container result(cleanup, 1, 100'000);
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 5u, 5u, 5 * (sequence_size - 2 + 1) });
}

TYPED_TEST(CountKmers, CountKmersAppendParallel) {
    using Container
            = common::SortedMultiset<TypeParam, uint8_t, Vector<std::pair<TypeParam, uint8_t>>>;
    Container result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    ASSERT_EQ(3u, result.data().size());

    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'B')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'B')), 2, &result,
            { 1 });
#if _DNA_GRAPH
    ASSERT_EQ(3u, result.data().size());
#else
    ASSERT_EQ(6u, result.data().size());
#endif
}

}  // namespace
