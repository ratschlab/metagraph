#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "dbg_succinct_merge.hpp"
#include "dbg_succinct_construct.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F,
                                 Config::StateType state) {
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << *graph->last;
    EXPECT_EQ(last, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    ostr << *(graph->W);
    EXPECT_EQ(W, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    for (size_t i = 0; i < graph->F.size(); ++i) {
        ostr << graph->F[i] << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state;
}


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F) {
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::STAT);
}


TEST(Construct, GraphDefaultConstructor) {
    DBG_succ *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ;
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ();
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        DBG_succ graph;
    });

    ASSERT_NO_THROW({
        DBG_succ graph();
    });
}

TEST(Construct, ConstructionEQAppendingSimplePath) {
    for (size_t k = 1; k < 80; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(100, 'A'));

        EXPECT_EQ(constructed, appended);
    }
}

TEST(Construct, ConstructionEQAppendingTwoPaths) {
    for (size_t k = 1; k < 80; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A'),
                                std::string(50, 'B') });
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        appended.add_sequence(std::string(100, 'A'));
        appended.add_sequence(std::string(50, 'B'));

        EXPECT_EQ(constructed, appended);
    }
}

TEST(Construct, ConstructionEQAppending) {
    for (size_t k = 1; k < 80; ++k) {
        std::vector<std::string> input_data = {
            "ACAGCTAGCTAGCTAGCTAGCTG",
            "ATATTATAAAAAATTTTAAAAAA",
            "ATATATTCTCTCTCTCTCATA",
            "GTGTGTGTGGGGGGCCCTTTTTTCATA",
        };
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads(input_data);
        DBG_succ constructed(&constructor);

        DBG_succ appended(k);
        for (const auto &sequence : input_data) {
            appended.add_sequence(sequence);
        }

        EXPECT_EQ(constructed, appended);
    }
}

TEST(DBGSuccinct, EmptyGraph) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ");
    delete graph;
}

TEST(DBGSuccinct, SwitchState) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::STAT);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    delete graph;
}

TEST(DBGSuccinct, AddSequenceFast) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    std::vector<std::string> names;

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
        names.emplace_back(read_stream->name.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    EXPECT_EQ(4llu, names.size());
    EXPECT_EQ("1", names[0]);
    EXPECT_EQ("2", names[1]);
    EXPECT_EQ("3", names[2]);
    EXPECT_EQ("4", names[3]);

    DBG_succ *graph = new DBG_succ(&constructor);

    //test graph construction
    test_graph(graph, "00011101101111111111111",
                      "00131124434010141720433",
                      "0 3 11 13 17 22 ");
    delete graph;
}

TEST(Construct, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    EXPECT_EQ(outgoing_edges[1], graph->outgoing(1, DBG_succ::encode('$')));

    for (size_t i = 1; i < graph->get_W().size(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBG_succ::encode('$')) {
            EXPECT_EQ(outgoing_edges[i], graph->outgoing(i, graph->get_W(i)))
                << "Edge index: " << i;

            EXPECT_EQ(
                graph->succ_last(i),
                graph->incoming(graph->outgoing(i, graph->get_W(i)),
                                graph->get_minus_k_value(i, graph->get_k() - 1).first)
            );
            for (TAlphabet c = 0; c < DBG_succ::alph_size; ++c) {
                uint64_t node_idx = graph->incoming(i, c);
                if (node_idx) {
                    EXPECT_EQ(
                        graph->succ_last(i),
                        graph->outgoing(node_idx, graph->get_node_last_value(i))
                    );
                }
            }
        }

        //test FM index property
        EXPECT_TRUE(graph->get_last(graph->fwd(i)));
        if (graph->get_W(i)) {
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_value(graph->fwd(i)));
        }
    }

    delete graph;
}

TEST(DBGSuccinct, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(&constructor);

    graph->serialize(test_dump_basename);

    DBG_succ loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";
    EXPECT_FALSE(DBG_succ() == loaded_graph);

    delete graph;
}

TEST(DBGSuccinct, AddSequenceSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A'));
        EXPECT_EQ(k + 1, graph.num_nodes());
        EXPECT_EQ(k + 2, graph.num_edges());
    }
}

TEST(DBGSuccinct, AddSequenceBugRevealingTestcase) {
    DBG_succ graph(1);
    graph.add_sequence("CTGAG", false);
}

TEST(DBGSuccinct, AddSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        graph.add_sequence("GAAC");
        EXPECT_EQ(11u, graph.num_nodes());
        EXPECT_EQ(14u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("AACG");
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC");
        graph.add_sequence("GACT");
        graph.add_sequence("ACTA");
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC", true);
        graph.add_sequence("AACG", true);
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(7u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(7u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAACT", true);
        graph.add_sequence("AAATG", true);
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequenceAnyKmerSize) {
    for (size_t k = 1; k < 10; ++k) {
        {
            DBG_succ graph(k);
            graph.add_sequence("AAAC", true);
            graph.add_sequence("AACG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", true);
            graph.add_sequence("AAATG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
        }
    }
}

void test_pred_kmer(const DBG_succ &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::deque<TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(), DBG_succ::encode);
    EXPECT_EQ(expected_idx, graph.pred_kmer(kmer)) << kmer_s << std::endl << graph;
}

TEST(DBGSuccinct, PredKmer) {
    {
        DBG_succ graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAAAA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("ACACA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(2);
        graph.add_sequence("ATAATATCC");
        graph.add_sequence("ATACGC");
        graph.add_sequence("ATACTC");
        graph.add_sequence("ATACTA");
        graph.add_sequence("CATT");
        graph.add_sequence("CCC");
        graph.add_sequence("GGGC");
        graph.add_sequence("GGTGTGAC");
        graph.add_sequence("GGTCT");
        graph.add_sequence("GGTA");

        test_pred_kmer(graph, "$$", 4);
        test_pred_kmer(graph, "$A", 5);
        test_pred_kmer(graph, "$T", 26);
        test_pred_kmer(graph, "AT", 29);
        test_pred_kmer(graph, "TT", 35);
        test_pred_kmer(graph, "NT", 35);
        test_pred_kmer(graph, "TN", 35);
    }
}

TEST(DBGSuccinct, PredKmerRandomTest) {
    for (size_t k = 1; k < 8; ++k) {
        DBG_succ graph(k);

        for (size_t p = 0; p < 10; ++p) {
            size_t length = rand() % 400;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, true);
        }

        auto all_kmer_str = utils::generate_strings("ACGT", k);
        for (size_t i = 1; i < k; ++i) {
            auto kmer_str_suffices = utils::generate_strings("ACGT", i);
            for (size_t j = 0; j < kmer_str_suffices.size(); ++j) {
                all_kmer_str.push_back(std::string(k - i, '$')
                                        + kmer_str_suffices[j]);
            }
        }

        for (const auto &kmer_str : all_kmer_str) {
            std::deque<TAlphabet> kmer(kmer_str.size());
            std::transform(kmer_str.begin(), kmer_str.end(),
                           kmer.begin(), DBG_succ::encode);

            uint64_t lower_bound = graph.pred_kmer(kmer);

            EXPECT_FALSE(
                utils::colexicographically_greater(
                    graph.get_node_seq(lower_bound), kmer
                )
            ) << graph
              << "kmer: " << kmer_str << std::endl
              << "lower bound: " << lower_bound << std::endl
              << "which is: " << graph.get_node_str(lower_bound) << std::endl;

            if (lower_bound < graph.get_W().size() - 1) {
                EXPECT_TRUE(
                    utils::colexicographically_greater(
                        graph.get_node_seq(lower_bound + 1), kmer
                    )
                );
            }
        }
    }
}

TEST(DBGSuccinct, FindSequence) {
    for (size_t k = 1; k < 10; ++k) {
        SequenceGraph *graph = new DBG_succ(k);

        graph->add_sequence(std::string(100, 'A'));

        uint64_t index = 0;
        graph->align(std::string(k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 2, index);
        graph->align(std::string(2 * k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 2, index);

        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 1));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 0));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 1));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 0));

        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 1));
        EXPECT_TRUE(graph->find(std::string(k, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 1));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'B'), 0));

        delete graph;
    }
}

TEST(DBGSuccinct, Traversals) {
    for (size_t k = 1; k < 10; ++k) {
        SequenceGraph *graph = new DBG_succ(k);

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        uint64_t it = 0;
        graph->align(std::string(k, 'A'), [&](uint64_t i) { it = i; });
        ASSERT_EQ(k + 3, it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it + 1, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it + 1, 'A'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse_back(it + 1, 'G'));

        delete graph;
    }
}


void sequence_to_kmers(std::vector<TAlphabet>&& seq,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix);

TEST(ExtractKmers, ExtractKmersEmptySuffix) {
    for (size_t k = 1; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length < k; ++length) {
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, {});
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 700; ++length) {
            result.clear();
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, {});
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers, ExtractKmersWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 1; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length < 500; ++length) {
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ 6 });
    for (size_t k = 1; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length < k; ++length) {
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_EQ(length - k, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }
    }
}

TEST(ExtractKmers, ExtractKmersWithFilteringTwo) {
    std::vector<TAlphabet> suffix = { 6, 1 };
    for (size_t k = 2; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length <= k; ++length) {
            sequence_to_kmers(std::vector<TAlphabet>(length, 6), k, &result, suffix);
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 500; ++length) {
            result.clear();

            std::vector<TAlphabet> sequence(length, 6);
            sequence[k - 1] = 1;

            sequence_to_kmers(std::move(sequence), k, &result, suffix);
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers, ExtractKmersAppend) {
    std::vector<KMer> result;

    sequence_to_kmers(std::vector<TAlphabet>(500, 6), 1, &result, {});
    ASSERT_EQ(499u, result.size());

    sequence_to_kmers(std::vector<TAlphabet>(500, 6), 1, &result, {});
    ASSERT_EQ(499u * 2, result.size());
}


void sequence_to_kmers(const std::string &sequence,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix);

TEST(ExtractKmers, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length < k; ++length) {
            sequence_to_kmers(std::string(length, 'N'), k, &result, {});
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            sequence_to_kmers(std::string(length, 'N'), k, &result, {});
            // NNN -> $NNN$
            ASSERT_EQ(length - k + 2, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TEST(ExtractKmers, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 2; k < 85; ++k) {
        std::vector<KMer> result;

        for (size_t length = 0; length < k; ++length) {
            sequence_to_kmers(std::string(length, 'N'), k, &result, { 5, 5 });
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();
            sequence_to_kmers(std::string(length, 'N'), k, &result, { 0, 0 });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            sequence_to_kmers(std::string(length, 'N'), k, &result, { 0, 5 });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            sequence_to_kmers(std::string(length, 'N'), k, &result, { 5, 5 });
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length <= k; ++length) {
            sequence_to_kmers(std::string(length, 'N'), k, &result, { 5, 1 });
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k + 1; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'N');
            sequence[k - 1] = 'A';

            sequence_to_kmers(sequence, k, &result, { 5, 1 });
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TEST(ExtractKmers, ExtractKmersFromStringAppend) {
    std::vector<KMer> result;

    sequence_to_kmers(std::string(500, 'A'), 1, &result, {});
    ASSERT_EQ(501u, result.size());

    sequence_to_kmers(std::string(500, 'A'), 1, &result, {});
    ASSERT_EQ(501u * 2, result.size());
}


typedef std::function<void(const std::string&)> CallbackRead;

void extract_kmers(std::function<void(CallbackRead)> generate_reads,
                   size_t k,
                   std::vector<KMer> *kmers,
                   size_t *end_sorted,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex *mutex,
                   bool remove_redundant = true,
                   size_t preallocated_size = 0);

void sequence_to_kmers_parallel_wrapper(std::vector<std::string> *reads,
                                        size_t k,
                                        std::vector<KMer> *kmers,
                                        const std::vector<TAlphabet> &suffix,
                                        std::mutex *mutex,
                                        bool remove_redundant) {
    size_t end_sorted = 0;
    kmers->reserve(100'000);
    extract_kmers([reads](CallbackRead callback) {
            for (auto &&read : *reads) {
                callback(std::move(read));
            }
        },
        k, kmers, &end_sorted, suffix,
        1, false, mutex, remove_redundant, 100'000
    );
    delete reads;
}

TEST(ExtractKmers, ExtractKmersAppendParallel) {
    std::vector<KMer> result;
    std::mutex mu;
    size_t sequence_size = 500;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, &mu, false
    );
    ASSERT_EQ((sequence_size + 1) * 5, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, &mu, false
    );
    ASSERT_EQ((sequence_size + 1) * 10, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'A')),
        1, &result, {}, &mu, false
    );
    ASSERT_EQ((sequence_size + 1) * 15, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, {}, &mu, false
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(sequence_size, 'B')),
        1, &result, { 1, }, &mu, false
    );
    ASSERT_EQ((sequence_size + 1) * 20, result.size());
}

TEST(ExtractKmers, ExtractKmersParallelRemoveRedundant) {
    std::vector<KMer> result;
    std::mutex mu;

    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        1, &result, {}, &mu, true
    );
    // $A, AA, A$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, {}, &mu, true
    );
    // $AA, AAA, AA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 0 }, &mu, true
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        2, &result, { 1 }, &mu, true
    );
    // $$A, $AA, AAA, AA$
    ASSERT_EQ(4u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, {}, &mu, true
    );
    // $AAA, AAAA, AAA$
    ASSERT_EQ(3u, result.size());

    result.clear();
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 0 }, &mu, true
    );
    sequence_to_kmers_parallel_wrapper(
        new std::vector<std::string>(5, std::string(500, 'A')),
        3, &result, { 1 }, &mu, true
    );
    // $$$A, $$AA, $AAA, AAAA, AAA$
    ASSERT_EQ(5u, result.size());
}
