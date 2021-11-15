#include <cstdio>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include <gtest/gtest.h>

#define protected public
#define private public

#include "common/seq_tools/reverse_complement.hpp"
#include "common/sorted_sets/sorted_set.hpp"
#include "common/sorted_sets/sorted_multiset.hpp"
#include "common/sorted_sets/sorted_multiset_disk.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "kmer/kmer_collector.hpp"
#include "tests/utils/gtest_patch.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph::boss;

using mtg::kmer::KmerExtractorBOSS;

const std::string test_data_dir = TEST_DATA_DIR;
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

#define kMaxK ( 256 / KmerExtractorBOSS::bits_per_char )


template <typename Kmer>
class CollectKmers : public ::testing::Test { };

template <typename Kmer>
class CountKmers : public ::testing::Test { };

typedef ::testing::Types<kmer::KMerBOSS<uint64_t, KmerExtractorBOSS::bits_per_char>,
                         kmer::KMerBOSS<sdsl::uint128_t, KmerExtractorBOSS::bits_per_char>,
                         kmer::KMerBOSS<sdsl::uint256_t, KmerExtractorBOSS::bits_per_char>> KmerTypes;

TYPED_TEST_SUITE(CollectKmers, KmerTypes);
TYPED_TEST_SUITE(CountKmers, KmerTypes);


TEST(BOSSConstruct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));

        for (bool weighted : { false, true }) {
            BOSSConstructor constructor(k, false, weighted ? 8 : 0);
            constructor.add_sequences({ std::string(100, 'A') });
            BOSS constructed(&constructor);

            EXPECT_EQ(constructed, appended);
        }
    }
}

TEST(BOSSConstruct, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        appended.add_sequence(std::string(100, 'A'));
        appended.add_sequence(std::string(50, 'B'));

        for (bool weighted : { false, true }) {
            BOSSConstructor constructor(k, false, weighted ? 8 : 0);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(50, 'B') });
            BOSS constructed(&constructor);

            EXPECT_EQ(constructed, appended);
        }
    }
}

TEST(BOSSConstruct, ConstructionLowerCase) {
    for (size_t k = 1; k < kMaxK; ++k) {
        for (bool weighted : { false, true }) {
            BOSSConstructor constructor_first(k, false, weighted ? 8 : 0);
            constructor_first.add_sequences({ std::string(100, 'A'),
                                              std::string(50, 'C') });
            BOSS first(&constructor_first);

            BOSSConstructor constructor_second(k, false, weighted ? 8 : 0);
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
}

TEST(BOSSConstruct, ConstructionDummySentinel) {
    for (size_t k = 1; k < kMaxK; ++k) {
        for (bool weighted : { false, true }) {
            BOSSConstructor constructor_first(k, false, weighted ? 8 : 0);
            constructor_first.add_sequences({ std::string(100, 'N'),
                                              std::string(50, '$') });
            BOSS first(&constructor_first);

            BOSSConstructor constructor_second(k, false, weighted ? 8 : 0);
            constructor_second.add_sequences({ std::string(100, 'N'),
                                               std::string(50, '.') });
            BOSS second(&constructor_second);

            EXPECT_TRUE(first.equals_internally(second));
        }
    }
}

TEST(BOSSConstruct, ConstructionEQAppending) {
    std::vector<std::string> input_data = {
        "ACAGCTAGCTAGCTAGCTAGCTG",
        "ATATTATAAAAAATTTTAAAAAA",
        "ATATATTCTCTCTCTCTCATA",
        "GTGTGTGTGGGGGGCCCTTTTTTCATA",
    };
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        for (const auto &sequence : input_data) {
            appended.add_sequence(sequence);
        }
        for (bool weighted : { false, true }) {
            for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
                BOSSConstructor constructor(k, false, weighted ? 8 : 0, "", 1,
                                            20000, container);
                constructor.add_sequences(std::vector<std::string>(input_data));
                BOSS constructed(&constructor);

                EXPECT_EQ(constructed, appended);
            }
        }
    }
}

TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeight) {
    std::vector<std::string> input_data = {
        "ACAGCTAGCTAGCTAGCTAGCTG",
        "ATATTATAAAAAATTTTAAAAAA",
        "ATATATTCTCTCTCTCTCATA",
        "GTGTGTGTGGGGGGCCCTTTTTTCATA",
    };
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            BOSSConstructor constructor(k, false, 8, "", 1, 20000, container);
            constructor.add_sequences(std::vector<std::string>(input_data));

            BOSS constructed;
            sdsl::int_vector_buffer<> weights;
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

TEST(WeightedBOSSConstruct, ConstructionDummyKmersZeroWeightChunks) {
    std::vector<std::string> input_data = {
        "ACAGCTAGCTAGCTAGCTAGCTG",
        "ATATTATAAAAAATTTTAAAAAA",
        "ATATATTCTCTCTCTCTCATA",
        "GTGTGTGTGGGGGGCCCTTTTTTCATA",
    };
    for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
        for (size_t k = 1; k < kMaxK; ++k) {
            BOSS constructed(k);

            auto constructor
                    = IBOSSChunkConstructor::initialize(k, false, 8, "", 1, 20000, container);

            for (auto &&sequence : input_data) {
                constructor->add_sequence(std::move(sequence));
            }

            BOSS::Chunk chunk = constructor->build_chunk();

            chunk.initialize_boss(&constructed);
            sdsl::int_vector_buffer<> weights = chunk.get_weights();

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

#if ! _PROTEIN_GRAPH
TEST(BOSSConstruct, ConstructionEQAppendingCanonical) {
    std::vector<std::string> input_data = {
        "ACAGCTAGCTAGCTAGCTAGCTG",
        "ATATTATAAAAAATTTTAAAAAA",
        "ATATATTCTCTCTCTCTCATA",
        "GTGTGTGTGGGGGGCCCTTTTTTCATA",
    };
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        for (auto &sequence : input_data) {
            appended.add_sequence(sequence);
            reverse_complement(sequence.begin(), sequence.end());
            appended.add_sequence(sequence);
        }
        for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
            for (bool weighted : { false, true }) {
                BOSSConstructor constructor(k, true, weighted ? 8 : 0, "", 1,
                                            20'000, container);
                constructor.add_sequences(std::vector<std::string>(input_data));
                BOSS constructed(&constructor);

                EXPECT_EQ(constructed, appended);
            }
        }
    }
}
#endif

TEST(BOSSConstruct, ConstructionLong) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        appended.add_sequence(std::string(k + 1, 'A'));

        for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
            for (bool weighted : { false, true }) {
                BOSSConstructor constructor(k, false, weighted ? 8 : 0, "", 1,
                                            20'000, container);
                constructor.add_sequences({ std::string(k + 1, 'A') });
                BOSS constructed(&constructor);

                EXPECT_EQ(constructed, appended);
                ASSERT_TRUE(constructed.num_nodes() > 1u);
            }
        }
    }
}

TEST(BOSSConstruct, ConstructionShort) {
    for (size_t k = 1; k < kMaxK; ++k) {
        BOSS appended(k);
        appended.add_sequence(std::string(k, 'A'));

        for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
            for (bool weighted : { false, true }) {
                BOSSConstructor constructor(k, false, weighted ? 8 : 0, "", 1,
                                            20'000, container);
                constructor.add_sequences({ std::string(k, 'A') });
                BOSS constructed(&constructor);

                EXPECT_EQ(constructed, appended);
                ASSERT_EQ(1u, constructed.num_nodes());
            }
        }
    }
}

TEST(BOSSConstruct, ConstructionFromChunks) {
    for (size_t k = 1; k < kMaxK; k += 6) {
        BOSS boss_dynamic(k);
        boss_dynamic.add_sequence(std::string(100, 'A'));
        boss_dynamic.add_sequence(std::string(100, 'C'));
        boss_dynamic.add_sequence(std::string(100, 'T') + "A"
                                        + std::string(100, 'G'));

        for (auto container : { kmer::ContainerType::VECTOR, kmer::ContainerType::VECTOR_DISK }) {
            for (size_t suffix_len = 0; suffix_len < std::min(k, (size_t)3u); ++suffix_len) {
                for (bool weighted : { false, true }) {
                    for (size_t num_threads : { 1, 4 }) {
                        BOSS::Chunk graph_data;

                        for (const std::string &suffix : KmerExtractorBOSS::generate_suffixes(suffix_len)) {
                            std::unique_ptr<IBOSSChunkConstructor> constructor(
                                    IBOSSChunkConstructor::initialize(k, false, weighted ? 8 : 0,
                                                                      suffix, num_threads, 20000, container));

                            constructor->add_sequence(std::string(100, 'A'));
                            constructor->add_sequence(std::string(100, 'C'));
                            constructor->add_sequence(std::string(100, 'T') + "A"
                                                            + std::string(100, 'G'));

                            BOSS::Chunk next_block = constructor->build_chunk();
                            if (graph_data.size()) {
                                graph_data.extend(next_block);
                            } else {
                                graph_data = std::move(next_block);
                            }
                        }

                        BOSS boss;
                        graph_data.initialize_boss(&boss);

                        EXPECT_EQ(boss_dynamic, boss);
                    }
                }
            }
        }
    }
}

template <typename KMER, class Container>
using Collector = typename mtg::kmer::KmerCollector<KMER, KmerExtractorBOSS, Container>;

// TODO: k is node length
template <typename KMER>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        common::SortedSet<typename KMER::WordType> *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    kmer::extract_kmers<KMER, KmerExtractorBOSS, common::SortedSet<typename KMER::WordType>>(
        [reads](kmer::CallString callback) {
            std::for_each(reads->begin(), reads->end(), callback);
        },
        k, Collector<KMER, common::SortedSet<typename KMER::WordType>>::BASIC, kmers, suffix
    );
    delete reads;
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallelReserved) {
    common::SortedSet<typename TypeParam::WordType> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, 100'000
    );
#if _DNA_GRAPH
    ASSERT_EQ(1u, result.data().size());
#else
    ASSERT_EQ(2u, result.data().size());
#endif

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, 100'000
    );
#if _DNA_GRAPH
    ASSERT_EQ(1u, result.data().size());
#else
    ASSERT_EQ(2u, result.data().size());
#endif
}

TYPED_TEST(CollectKmers, CollectKmersAppendParallel) {
    common::SortedSet<typename TypeParam::WordType> result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 0
    );
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper<TypeParam>(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, 0
    );
#if _DNA_GRAPH
    ASSERT_EQ(1u, result.data().size());
#else
    ASSERT_EQ(2u, result.data().size());
#endif
}

// TODO: k is node length
template <typename KMER, typename Container>
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        Container *kmers,
                                        const std::vector<KmerExtractorBOSS::TAlphabet> &suffix) {
    static_assert(std::is_same_v<typename Container::key_type, typename KMER::WordType>);
    // kmers->try_reserve(reserved_capacity);
    kmer::count_kmers<KMER, KmerExtractorBOSS, Container>(
        [reads](kmer::CallStringCount callback) {
            for (const auto &read : *reads) {
                callback(read, 1);
            }
        },
        k, Collector<KMER, Container>::BASIC, kmers, suffix
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

template <typename KMER, typename Container>
void check_counts() {
    Container result(1, 100'000);
    size_t sequence_size = 500;
    size_t max_value = std::numeric_limits<typename Container::count_type>::max();
    const size_t five_times = std::min(5 * (sequence_size - 2 + 1), max_value);
    const size_t ten_times = std::min(10 * (sequence_size - 2 + 1), max_value);
    const size_t fifteen_times = std::min(15 * (sequence_size - 2 + 1), max_value);

    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { five_times }); // AA - #five_times


    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { ten_times }); // AA - #ten_times


    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { fifteen_times }); // AA - #fifteen_times


    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result, {});
    // AA - fifteen_times, CC - five_times
    assert_contents(result, { fifteen_times, five_times });


    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result, { 1 });
    // same as before, as no k-mers end with 'A'
    assert_contents(result, { fifteen_times, five_times });


    sequence_to_kmers_parallel_wrapper<KMER, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'C')), 2, &result, { 0 });
    // $C - 5 (matches the given suffix), AA - fifteen_times, CC - five_times
    assert_contents(result, { 5u, fifteen_times, five_times });
}

TYPED_TEST(CountKmers, CountKmers8bits) {
    using Container = common::SortedMultiset<typename TypeParam::WordType, uint8_t>;
    check_counts<TypeParam, Container>();
}

TYPED_TEST(CountKmers, CountKmers32bits) {
    using Container = common::SortedMultiset<typename TypeParam::WordType, uint32_t>;
    check_counts<TypeParam, Container>();
}

TYPED_TEST(CountKmers, CountKmers8bitsDisk) {
    using Container = common::SortedMultisetDisk<typename TypeParam::WordType, uint8_t>;
    Container result(1, 100'000);
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 255u });
}

TYPED_TEST(CountKmers, CountKmers32bitsDisk) {
    using Container = common::SortedMultisetDisk<typename TypeParam::WordType, uint32_t>;
    Container result(1, 100'000);
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    assert_contents(result, { 5 * (sequence_size - 2 + 1) });
}

TYPED_TEST(CountKmers, CountKmersAppendParallel) {
    using Container = common::SortedMultiset<typename TypeParam::WordType, uint8_t>;
    Container result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'A')), 2, &result, {});
    // AA
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'B')), 2, &result, {});
    sequence_to_kmers_parallel_wrapper<TypeParam, Container>(
            new std::vector<std::string>(5, std::string(sequence_size, 'B')), 2, &result, { 1 });
#if _DNA_GRAPH
    ASSERT_EQ(1u, result.data().size());
#else
    // AA, BB
    ASSERT_EQ(2u, result.data().size());
#endif
}

}  // namespace
