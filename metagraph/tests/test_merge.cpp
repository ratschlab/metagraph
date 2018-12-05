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
#include "dbg_construct.hpp"
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
                sequence[s] = random.alphabet[1 + rand() % 4];
            }
            random.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = random.alphabet[1 + rand() % 4];
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

TEST(DBGSuccinctMerge, TraversalMergeDisconnectedGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        DBGSuccConstructor constructor_first(k);
        constructor_first.add_sequences({ std::string(100, 'A') });
        DBG_succ first(&constructor_first);
        ASSERT_EQ(2u, first.num_edges());

        DBGSuccConstructor constructor_second(k);
        constructor_second.add_sequences({ std::string(50, 'C'),
                                           std::string(60, 'G') });
        DBG_succ second(&constructor_second);
        ASSERT_EQ(3u, second.num_edges());

        DBGSuccConstructor constructor_third(k);
        constructor_third.add_sequences({ std::string(100, 'A'),
                                          std::string(50, 'C'),
                                          std::string(60, 'G') });
        DBG_succ result(&constructor_third);
        ASSERT_EQ(4u, result.num_edges());

        first.switch_state(Config::DYN);
        first.merge(second);

        EXPECT_EQ(result, first);
    }
}

TEST(DBGSuccinctMerge, ParallelMergeEmptyGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k);
        DBG_succ second(k);

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);

        std::unique_ptr<DBG_succ::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);

        std::unique_ptr<DBG_succ::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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
        DBG_succ first(k);
        DBG_succ second(k);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));

        std::vector<const DBG_succ*> graphs = { &first, &second };

        DBG_succ *merged = merge::merge(graphs);

        std::unique_ptr<DBG_succ::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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

        std::unique_ptr<DBG_succ::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, 1, 1)
        };
        chunk->serialize(test_data_dir + "/1");
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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

        {
            std::unique_ptr<DBG_succ::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 0, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/1");
        }
        {
            std::unique_ptr<DBG_succ::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 1, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/2");
        }
        {
            std::unique_ptr<DBG_succ::Chunk> chunk {
                merge::merge_blocks_to_chunk(graphs, 2, 3, 1, 1)
            };
            chunk->serialize(test_data_dir + "/3");
        }
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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

        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
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

TEST(DBGSDMerge, DumpedChunked) {
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

        ASSERT_TRUE(full.equals_internally(*chunked, true)) << k;
    }
}

TEST(DBGSDMerge, DumpedChunkedCanonical) {
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

        ASSERT_TRUE(full.equals_internally(*chunked, true)) << k;
    }
}

TEST(DBGSDMerge, ParallelDumpedChunked) {
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

        ASSERT_TRUE(full.equals_internally(*chunked, true)) << k;
    }
}

TEST(DBGSDMerge, ParallelDumpedChunkedCanonical) {
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

        ASSERT_TRUE(full.equals_internally(*chunked, true)) << k;
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

        DBG_succ *merged = merge::merge(graphs);

        std::unique_ptr<DBG_succ::Chunk> chunk {
            merge::merge_blocks_to_chunk(graphs, 0, 1, num_threads, num_bins_per_thread)
        };
        chunk->serialize(test_data_dir + "/1");
        DBG_succ *chunked_merged = DBG_succ::Chunk::build_graph_from_chunks(
            { test_data_dir + "/1" }
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
