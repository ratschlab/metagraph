#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif
#endif

#define protected public
#define private public

#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

typedef KMer<uint64_t, KmerExtractor2Bit::kLogSigma> KMER;
const int kMaxK = sizeof(KMER) * 8 / KmerExtractor2Bit::kLogSigma - 1;

const KmerExtractor2Bit kmer_extractor;


TEST(Construct_SD_64, ConstructionNEAppendingSimplePath) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGSD constructed(&constructor);
        EXPECT_EQ(1u, constructed.num_nodes());

        DBGSD appended(k);
        ASSERT_DEATH(appended.add_sequence(std::string(100, 'A')), "");

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
        ASSERT_DEATH(appended.add_sequence(std::string(100, 'A')), "");
        ASSERT_DEATH(appended.add_sequence(std::string(50, 'C')), "");

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
            kmer_extractor.sequence_to_kmers(str, k, {}, &kmers);
        }
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        EXPECT_EQ(kmers.size(), constructed.num_nodes());

        DBGSD appended(k);
        EXPECT_NE(constructed, appended);
    }
}

TEST(Construct_SD_64, ConstructionEQAppendingCanonical) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        DBGSDConstructor constructor(k, true);
        constructor.add_sequences(input_data);
        DBGSD constructed(&constructor);

        Vector<KmerExtractor2Bit::Kmer64> kmers;
        for (const auto &str : input_data) {
            kmer_extractor.sequence_to_kmers(str, k, {}, &kmers);
            auto rev_str = str;
            reverse_complement(rev_str.begin(), rev_str.end());
            kmer_extractor.sequence_to_kmers(rev_str, k, {}, &kmers);
        }

        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        EXPECT_EQ(kmers.size(), constructed.num_nodes());

        DBGSD appended(k);
        EXPECT_NE(constructed, appended);
    }
}

using TAlphabet = KmerExtractor2Bit::TAlphabet;

TEST(ExtractKmersPacked_64, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        // NNN -> $NNN$
        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
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
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ kmer_extractor.encode('T') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<KMER> result;

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
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
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TT"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, { 0, kmer_extractor.encode('T') }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TT"), &result
            );

            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'T');
            sequence[k - 1] = 'A';

            kmer_extractor.sequence_to_kmers(
                sequence, k, kmer_extractor.encode("TA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmersPacked_64, ExtractKmersFromStringAppend) {
    Vector<KMER> result;

    kmer_extractor.sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u, result.size());

    kmer_extractor.sequence_to_kmers(
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
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    ASSERT_EQ(2u, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, mu_resize, mu, false, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, mu_resize, mu, false, 0
    );
    sort_and_remove_duplicates(&result, 1);
    // #if _DNA4_GRAPH
        // B->A in 2Bit mode
        ASSERT_EQ(2u, result.size());
    // #else
    //     ASSERT_EQ(3u, result.size());
    // #endif
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

TEST(DBGSDMergeChunks, Chunked) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "A"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "C"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "G"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "T"));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::unique_ptr<DBGSD::Chunk>> chunks;

        for (size_t i = 0; i < constructors.size(); ++i) {
            chunks.emplace_back(constructors[i]->build_chunk());
            ASSERT_TRUE(chunks.back().get());
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(chunks)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, ChunkedCanonical) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "A"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "C"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "G"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "T"));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::unique_ptr<DBGSD::Chunk>> chunks;

        for (size_t i = 0; i < constructors.size(); ++i) {
            chunks.emplace_back(constructors[i]->build_chunk());
            ASSERT_TRUE(chunks.back().get());
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(chunks, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, ParallelChunked) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "A", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "C", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "G", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "T", num_threads));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::unique_ptr<DBGSD::Chunk>> chunks;
        uint64_t chunk_size = 0;

        for (size_t i = 0; i < constructors.size(); ++i) {
            chunks.emplace_back(constructors[i]->build_chunk());
            ASSERT_TRUE(chunks.back().get());
            chunk_size += chunks.back()->num_set_bits();
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(chunks)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, ParallelChunkedCanonical) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "A", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "C", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "G", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "T", num_threads));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::unique_ptr<DBGSD::Chunk>> chunks;
        uint64_t chunk_size = 0;

        for (size_t i = 0; i < constructors.size(); ++i) {
            chunks.emplace_back(constructors[i]->build_chunk());
            ASSERT_TRUE(chunks.back().get());
            chunk_size += chunks.back()->num_set_bits();
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(chunks, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, DumpedChunked) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "A"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "C"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "G"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "T"));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::string> files;

        for (size_t i = 0; i < constructors.size(); ++i) {
            auto chunk = constructors[i]->build_chunk();
            ASSERT_TRUE(chunk);
            files.push_back(test_data_dir + "/chunks_to_merge"
                              + "." + std::to_string(i)
                              + "_" + std::to_string(4)
                              + ".dbgsdchunk");
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(files)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, DumpedChunkedCanonical) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "A"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "C"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "G"));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "T"));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::string> files;

        for (size_t i = 0; i < constructors.size(); ++i) {
            auto chunk = constructors[i]->build_chunk();
            ASSERT_TRUE(chunk);
            files.push_back(test_data_dir + "/chunks_to_merge"
                              + "." + std::to_string(i)
                              + "_" + std::to_string(4)
                              + ".dbgsdchunk");
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(files, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, ParallelDumpedChunked) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "A", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "C", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "G", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, false, "T", num_threads));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::string> files;
        uint64_t chunk_size = 0;

        for (size_t i = 0; i < constructors.size(); ++i) {
            auto chunk = constructors[i]->build_chunk();
            ASSERT_TRUE(chunk);
            chunk_size += chunk->num_set_bits();
            files.push_back(test_data_dir + "/chunks_to_merge"
                              + "." + std::to_string(i)
                              + "_" + std::to_string(4)
                              + ".dbgsdchunk");
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(files)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGSDMergeChunks, ParallelDumpedChunkedCanonical) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<ISDChunkConstructor>> constructors;
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "A", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "C", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "G", num_threads));
        constructors.emplace_back(ISDChunkConstructor::initialize(k, true, "T", num_threads));

        for (auto &constructor : constructors) {
            constructor->add_sequence("AAACT");
            constructor->add_sequence("ACTATG");
            constructor->add_sequence(std::string(50, 'C'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'A'));
            constructor->add_sequence(std::string(60, 'T'));
            constructor->add_sequence(std::string(60, 'G'));
            constructor->add_sequence(std::string(60, 'C'));
        }

        std::vector<std::string> files;
        uint64_t chunk_size = 0;

        for (size_t i = 0; i < constructors.size(); ++i) {
            auto chunk = constructors[i]->build_chunk();
            ASSERT_TRUE(chunk);
            chunk_size += chunk->num_set_bits();
            files.push_back(test_data_dir + "/chunks_to_merge"
                              + "." + std::to_string(i)
                              + "_" + std::to_string(4)
                              + ".dbgsdchunk");
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGSD> chunked{
            DBGSDConstructor::build_graph_from_chunks(files, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGSDConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGSD full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}
