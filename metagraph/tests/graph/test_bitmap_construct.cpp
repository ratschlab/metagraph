#include <stdio.h>
#include <string>
#include <sstream>
#include <filesystem>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"
#include "../test_helpers.hpp"


#define protected public
#define private public

#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"
#include "common/utils/string_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/sorted_sets/sorted_set.hpp"
#include "kmer/kmer_collector.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using mtg::kmer::KmerExtractor2Bit;

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

typedef kmer::KMer<uint64_t, KmerExtractor2Bit::bits_per_char> KMER;
const int kMaxK = sizeof(KMER) * 8 / KmerExtractor2Bit::bits_per_char - 1;

const KmerExtractor2Bit kmer_extractor;


TEST(DBGBitmapConstruct, ConstructionNEAppendingSimplePath) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGBitmap constructed(&constructor);
        EXPECT_EQ(1u, constructed.num_nodes());

        DBGBitmap appended(k);
        ASSERT_THROW(appended.add_sequence(std::string(100, 'A')), std::runtime_error);

        ASSERT_NE(constructed, appended);
    }
}

TEST(DBGBitmapConstruct, ConstructionNEAppendingTwoPaths) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(50, 'C') });
        DBGBitmap constructed(&constructor);
        EXPECT_EQ(2u, constructed.num_nodes());

        DBGBitmap appended(k);
        ASSERT_THROW(appended.add_sequence(std::string(100, 'A')), std::runtime_error);
        ASSERT_THROW(appended.add_sequence(std::string(50, 'C')), std::runtime_error);

        ASSERT_NE(constructed, appended);
    }
}

// TODO
/*
TEST(DBGBitmapConstruct, ConstructionLowerCase) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        DBGBitmapConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C') });
        DBGBitmap first(&constructor_first);

        DBGBitmapConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(100, 'a'),
                                           std::string(50, 'c') });
        DBGBitmap second(&constructor_second);

#if _DNA_CASE_SENSITIVE_GRAPH
        EXPECT_FALSE(first.equals_internally(second));
#else
        EXPECT_TRUE(first.equals_internally(second));
#endif
    }
}
*/

TEST(DBGBitmapConstruct, ConstructionEQAppending) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences(std::vector<std::string>(input_data));
        DBGBitmap constructed(&constructor);

        Vector<KmerExtractor2Bit::Kmer64> kmers;
        for (const auto &str : input_data) {
            kmer_extractor.sequence_to_kmers(str, k, {}, &kmers);
        }
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        EXPECT_EQ(kmers.size(), constructed.num_nodes());

        DBGBitmap appended(k);
        EXPECT_NE(constructed, appended);
    }
}

TEST(DBGBitmapConstruct, ConstructionEQAppendingCanonical) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        DBGBitmapConstructor constructor(k, true);
        constructor.add_sequences(std::vector<std::string>(input_data));
        DBGBitmap constructed(&constructor);

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

        DBGBitmap appended(k);
        EXPECT_NE(constructed, appended);
    }
}

TEST(DBGBitmapConstruct, ConstructionFromChunks) {
    for (bool canonical : { true, false }) {
        for (size_t k = 2; k <= kMaxK; k += 6) {
            std::vector<std::string> input_data = {
                "ACAGCTAGCTAGCTAGCTAGCTG",
                "ATATTATAAAAAATTTTAAAAAA",
                "ATATATTCTCTCTCTCTCATA",
                "GTGTGTGTGGGGGGCCCTTTTTTCATA",
                std::string(100, 'A'),
                std::string(100, 'C')
            };
            auto constructor = std::make_unique<DBGBitmapConstructor>(k, canonical);
            constructor->add_sequences(std::vector<std::string>(input_data));
            DBGBitmap reference(constructor.get());

            for (size_t suffix_len = 0; suffix_len <= k && suffix_len <= 4u; ++suffix_len) {
                std::vector<std::string> chunk_filenames;

                //one pass per suffix
                for (const std::string &suffix : KmerExtractor2Bit().generate_suffixes(suffix_len)) {
                    constructor.reset(new DBGBitmapConstructor(k, canonical, 0, suffix));

                    for (const auto &seq : input_data) {
                        constructor->add_sequence(std::string(seq));
                    }

                    std::unique_ptr<DBGBitmap::Chunk> chunk { constructor->build_chunk() };

                    chunk_filenames.push_back(
                        utils::join_strings({ test_data_dir + "/chunks_to_merge", suffix }, ".")
                            + DBGBitmap::kChunkFileExtension
                    );

                    std::filesystem::remove(chunk_filenames.back());

                    std::ofstream out(chunk_filenames.back(), std::ios::binary);
                    chunk->serialize(out);
                }

                ASSERT_TRUE(chunk_filenames.size());

                std::unique_ptr<DBGBitmap> graph(constructor->build_graph_from_chunks(chunk_filenames, canonical));

                for (const auto &filename : chunk_filenames) {
                    std::filesystem::remove(filename);
                }

                EXPECT_EQ(reference, *graph);
            }
        }
    }
}


// TODO: k is node length
void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        common::SortedSet<KMER::WordType> *kmers,
                                        const std::vector<KmerExtractor2Bit::TAlphabet> &suffix,
                                        size_t reserved_capacity) {
    kmers->try_reserve(reserved_capacity);
    using KMER_INT = typename KMER::WordType ;
    kmer::extract_kmers<KMER, KmerExtractor2Bit, common::SortedSet<KMER_INT>>(
        [reads](kmer::CallString callback) {
            std::for_each(reads->begin(), reads->end(), callback);
        },
        k, false, kmers, suffix
    );
    delete reads;
}

TEST(CollectKmers2Bit, ExtractKmersAppendParallelReserved) {
    common::SortedSet<typename KMER::WordType > result;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    EXPECT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    EXPECT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        2, &result, {}, 100'000
    );
    EXPECT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, 100'000
    );
    EXPECT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, 100'000
    );
    EXPECT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, 100'000
    );
    EXPECT_EQ(2u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, 100'000
    );
    EXPECT_EQ(2u, result.data().size());
}

TEST(CollectKmers2Bit, ExtractKmersAppendParallel) {
    common::SortedSet<typename KMER::WordType> result;
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
    ASSERT_EQ(1u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'C')),
        2, &result, { 1, }, 0
    );
    ASSERT_EQ(2u, result.data().size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, {}, 0
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        2, &result, { 1, }, 0
    );
    // #if _DNA4_GRAPH
        // B->A in 2Bit mode
        ASSERT_EQ(2u, result.data().size());
    // #else
    //     ASSERT_EQ(3u, result.data().size());
    // #endif
}

TEST(DBGBitmapMergeChunks, DumpedChunked) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<IBitmapChunkConstructor>> constructors;
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "A"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "C"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "G"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "T"));

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
            std::filesystem::remove(files.back());
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGBitmap> chunked{
            DBGBitmapConstructor::build_graph_from_chunks(files)
        };

        ASSERT_TRUE(chunked.get());

        DBGBitmapConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGBitmap full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGBitmapMergeChunks, DumpedChunkedCanonical) {
    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<IBitmapChunkConstructor>> constructors;
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "A"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "C"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "G"));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "T"));

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
            std::filesystem::remove(files.back());
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGBitmap> chunked{
            DBGBitmapConstructor::build_graph_from_chunks(files, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGBitmapConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGBitmap full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGBitmapMergeChunks, ParallelDumpedChunked) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<IBitmapChunkConstructor>> constructors;
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "A", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "C", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "G", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, false, 0, "T", num_threads));

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
            std::filesystem::remove(files.back());
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGBitmap> chunked{
            DBGBitmapConstructor::build_graph_from_chunks(files)
        };

        ASSERT_TRUE(chunked.get());

        DBGBitmapConstructor full_constructor(k);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGBitmap full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

TEST(DBGBitmapMergeChunks, ParallelDumpedChunkedCanonical) {
    const size_t num_threads = 4;

    for (size_t k = 2; k < 11; ++k) {
        std::vector<std::unique_ptr<IBitmapChunkConstructor>> constructors;
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "A", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "C", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "G", num_threads));
        constructors.emplace_back(IBitmapChunkConstructor::initialize(k, true, 0, "T", num_threads));

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
            std::filesystem::remove(files.back());
            std::ofstream file(files.back(), std::ios::binary);
            chunk->serialize(file);
            delete chunk;
        }

        std::unique_ptr<DBGBitmap> chunked{
            DBGBitmapConstructor::build_graph_from_chunks(files, true)
        };

        ASSERT_TRUE(chunked.get());

        DBGBitmapConstructor full_constructor(k, true);
        full_constructor.add_sequence("AAACT");
        full_constructor.add_sequence("ACTATG");
        full_constructor.add_sequence(std::string(50, 'C'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'A'));
        full_constructor.add_sequence(std::string(60, 'T'));
        full_constructor.add_sequence(std::string(60, 'G'));
        full_constructor.add_sequence(std::string(60, 'C'));

        DBGBitmap full(2);
        full_constructor.build_graph(&full);
        ASSERT_EQ(full.num_nodes(), chunk_size);
        ASSERT_EQ(full_constructor.constructor_->get_k(), full.get_k());

        ASSERT_TRUE(full.equals(*chunked, true)) << k;
    }
}

} // namespace
