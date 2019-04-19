#include <stdio.h>
#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "boss.hpp"
#include "boss_merge.hpp"
#include "boss_construct.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";


TEST(DBGSuccinctMerge, TraversalMergeWithEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        second.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEmptyRandomTest) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS random(k);

        for (size_t p = 0; p < 5; ++p) {
            size_t length = rand() % 10;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = random.alphabet[1 + rand() % 4];
            }
            random.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = random.alphabet[1 + rand() % 4];
            }
            random.add_sequence(sequence, true);
        }

        BOSS result(k);
        result.merge(random);

        EXPECT_EQ(random, result);
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEqualPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        first.merge(second);
        EXPECT_EQ(2 * k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(2 * k + 3, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeSinglePathWithTwo) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        first.merge(second);
        EXPECT_EQ(3 * k + 1, first.num_nodes());
        EXPECT_EQ(2 * k + 1, second.num_nodes());
        EXPECT_EQ(3 * k + 4, first.num_edges());
        EXPECT_EQ(2 * k + 3, second.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeTwoGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        first.add_sequence(std::string(50, 'C'));
        first.add_sequence(std::string(60, 'G'));
        first.add_sequence("AAAGT");
        second.add_sequence("AAACT", true);
        second.add_sequence("AAATG", true);
        second.add_sequence("ACTGA", true);
        BOSS merged(k);
        merged.merge(second);
        merged.merge(first);
        first.merge(second);
        EXPECT_EQ(first, merged);
    }
}

TEST(DBGSuccinctMerge, TraversalMergeDisconnectedGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'A') });
        BOSS first(&constructor_first);
        ASSERT_EQ(2u, first.num_edges());

        BOSSConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(50, 'C'),
                                           std::string(60, 'G') });
        BOSS second(&constructor_second);
        ASSERT_EQ(3u, second.num_edges());

        BOSSConstructor constructor_third(k);
        constructor_third.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C'),
                                          std::string(60, 'G') });
        BOSS result(&constructor_third);
        ASSERT_EQ(4u, result.num_edges());

        first.switch_state(Config::DYN);
        first.merge(second);

        EXPECT_EQ(result, first);
    }
}

TEST(DBGSuccinctMerge, ParallelMergeEmptyGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);

        std::vector<const BOSS*> graphs = { &first, &second };

        BOSS *merged = merge::merge(graphs);

        std::unique_ptr<BOSS::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1" }
        );

        first.merge(second);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelMergeTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));

        std::vector<const BOSS*> graphs = { &first, &second };

        BOSS *merged = merge::merge(graphs);

        std::unique_ptr<BOSS::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1" }
        );

        first.merge(second);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelMergeSinglePathWithTwo) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));

        std::vector<const BOSS*> graphs = { &first, &second };

        BOSS *merged = merge::merge(graphs);

        std::unique_ptr<BOSS::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1" }
        );

        first.merge(second);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelMergeThreeGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        BOSS third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const BOSS*> graphs = { &first, &second, &third };

        BOSS *merged = merge::merge(graphs);

        std::unique_ptr<BOSS::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1" }
        );

        first.merge(second);
        first.merge(third);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelChunkedMergeThreeGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        BOSS third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const BOSS*> graphs = { &first, &second, &third };

        BOSS *merged = merge::merge(graphs);

        {
            std::unique_ptr<BOSS::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 0, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/1");
        }
        {
            std::unique_ptr<BOSS::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 1, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/2");
        }
        {
            std::unique_ptr<BOSS::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 2, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/3");
        }
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1",
              test_data_dir + "/2",
              test_data_dir + "/3", }
        );

        first.merge(second);
        first.merge(third);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelDumpedChunkedMergeThreeGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS first(k);
        BOSS second(k);
        BOSS third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const BOSS*> graphs = { &first, &second, &third };

        BOSS *merged = merge::merge(graphs);
        size_t num_chunks = 3;

        std::vector<std::string> files;
        for (size_t i = 0; i < num_chunks; ++i) {
            auto chunk = merge::merge_blocks_to_chunk(graphs, i, num_chunks, 1, 1);
            ASSERT_TRUE(chunk);
            files.push_back(test_data_dir + "/chunks_to_merge"
                              + "." + std::to_string(i)
                              + "_" + std::to_string(num_chunks));
            chunk->serialize(files.back());
            delete chunk;
        }

        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            files
        );
        ASSERT_TRUE(chunked_merged);

        first.merge(second);
        first.merge(third);

        EXPECT_EQ(first, *merged);
        EXPECT_EQ(first, *chunked_merged);
        delete merged;
        delete chunked_merged;
    }
}

void random_testing_parallel_merge(size_t num_graphs, size_t num_sequences, size_t max_length,
                                   size_t num_threads, size_t num_bins_per_thread) {
    for (size_t k = 1; k < 10; ++k) {
        std::vector<const BOSS*> graphs(num_graphs, NULL);

        for (size_t i = 0; i < graphs.size(); ++i) {
            BOSS *component = new BOSS(k);

            for (size_t p = 0; p < num_sequences; ++p) {
                size_t length = rand() % max_length;
                std::string sequence(length, 'A');

                for (size_t s = 0; s < sequence.size(); ++s) {
                    sequence[s] = component->alphabet[1 + rand() % 4];
                }
                component->add_sequence(sequence, false);

                for (size_t s = 0; s < sequence.size(); ++s) {
                    sequence[s] = component->alphabet[1 + rand() % 4];
                }
                component->add_sequence(sequence, true);
            }
            graphs[i] = component;
        }

        BOSS *merged = merge::merge(graphs);

        std::unique_ptr<BOSS::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, num_threads, num_bins_per_thread)
        };
        chunk->serialize(test_data_dir + "/1");
        BOSS *chunked_merged = BOSS::Chunk::build_boss_from_chunks(
            { test_data_dir + "/1" }
        );

        BOSS result(k);
        for (size_t i = 0; i < graphs.size(); ++i) {
            result.merge(*graphs[i]);
        }

        ASSERT_EQ(result, *merged) << "The first merged graph is:\n"
                                   << *graphs[0];

        ASSERT_EQ(result, *chunked_merged) << "The first merged graph is:\n"
                                           << *graphs[0];

        for (size_t i = 0; i < graphs.size(); ++i) {
            delete graphs[i];
        }

        delete merged;
        delete chunked_merged;
    }
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_1_5_10_1_1) {
    random_testing_parallel_merge(1, 5, 10, 1, 1);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_1_5_10_1_3) {
    random_testing_parallel_merge(1, 5, 10, 1, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_1_5_4_1_300) {
    random_testing_parallel_merge(1, 5, 4, 1, 300);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_3_5_4_1_3) {
    random_testing_parallel_merge(3, 5, 4, 1, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_3_3_1_1_1) {
    random_testing_parallel_merge(3, 3, 1, 1, 1);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_5_5_10_1_3) {
    random_testing_parallel_merge(5, 5, 10, 1, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_15_5_4_1_3) {
    random_testing_parallel_merge(15, 5, 4, 1, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_15_5_4_20_3) {
    random_testing_parallel_merge(15, 5, 4, 20, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_15_10_30_20_3) {
    random_testing_parallel_merge(15, 10, 30, 20, 3);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_15_5_4_40_10) {
    random_testing_parallel_merge(15, 5, 4, 40, 10);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_15_5_20_39_10) {
    random_testing_parallel_merge(15, 5, 20, 39, 10);
}

TEST(DBGSuccinctMerge, ParallelMergeGraphsRandom_20_10_10_40_9) {
    random_testing_parallel_merge(20, 10, 10, 40, 9);
}
