#include <stdio.h>
#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "dbg_succinct_merge.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";


TEST(DBGSuccinctMerge, TraversalMergeWithEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k);
        DBG_succ second(k);
        second.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEmptyRandomTest) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ random(k);

        for (size_t p = 0; p < 5; ++p) {
            size_t length = rand() % 10;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            random.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            random.add_sequence(sequence, true);
        }

        DBG_succ result(k);
        result.merge(random);

        EXPECT_EQ(random, result);
    }
}

TEST(DBGSuccinctMerge, TraversalMergeEqualPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k);
        DBG_succ second(k);
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
        DBG_succ first(k);
        DBG_succ second(k);
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
        DBG_succ first(k);
        DBG_succ second(k);
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
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        first.add_sequence(std::string(50, 'C'));
        first.add_sequence(std::string(60, 'G'));
        first.add_sequence("AAAGT");
        second.add_sequence("AAACT", true);
        second.add_sequence("AAATG", true);
        second.add_sequence("ACTGA", true);
        DBG_succ merged(k);
        merged.merge(second);
        merged.merge(first);
        first.merge(second);
        EXPECT_EQ(first, merged);
    }
}

TEST(DBGSuccinctMerge, ParallelMergeEmptyGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k);
        DBG_succ second(k);

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1) }
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
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1) }
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
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1) }
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
        DBG_succ first(k);
        DBG_succ second(k);
        DBG_succ third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const DBG_succ*> graphs = { &first, &second, &third };

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1) }
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
        DBG_succ first(k);
        DBG_succ second(k);
        DBG_succ third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const DBG_succ*> graphs = { &first, &second, &third };

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 3, 1, 1),
              merge::merge_blocks_to_chunk(graphs, 1, 3, 1, 1), 
              merge::merge_blocks_to_chunk(graphs, 2, 3, 1, 1) }
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
        DBG_succ first(k);
        DBG_succ second(k);
        DBG_succ third(k);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const DBG_succ*> graphs = { &first, &second, &third };

        DBG_succ *merged = merge::merge(graphs);
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

        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            std::vector<DBG_succ::Chunk*>(num_chunks, NULL),
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
        std::vector<const DBG_succ*> graphs(num_graphs, NULL);

        for (size_t i = 0; i < graphs.size(); ++i) {
            DBG_succ *component = new DBG_succ(k);

            for (size_t p = 0; p < num_sequences; ++p) {
                size_t length = rand() % max_length;
                std::string sequence(length, 'A');

                for (size_t s = 0; s < sequence.size(); ++s) {
                    sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
                }
                component->add_sequence(sequence, false);

                for (size_t s = 0; s < sequence.size(); ++s) {
                    sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
                }
                component->add_sequence(sequence, true);
            }
            graphs[i] = component;
        }

        DBG_succ *merged = merge::merge(graphs);
        DBG_succ *chunked_merged = merge::build_graph_from_chunks(
            k,
            { merge::merge_blocks_to_chunk(graphs, 0, 1, num_threads, num_bins_per_thread) }
        );
        DBG_succ result(k);
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
