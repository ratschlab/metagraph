#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_hash.hpp"


TEST(DBGHash, SimplePathReconstruct) {
    DBGHash graph(4);
    graph.add_sequence("ACGTGCTA");

    for (size_t i = graph.first_edge(); i <= graph.last_edge(); ++i) {
        auto kmer = graph.get_node_kmer(i);
        auto label = graph.get_edge_label(i);

        if (!graph.is_dummy_edge(kmer + label)) {
            ASSERT_TRUE(graph.has_the_only_outgoing_edge(i)) << kmer + label;
            ASSERT_TRUE(graph.has_the_only_incoming_edge(i)) << kmer + label;
        }
    }
}

TEST(DBGHash, SimplePathOutDegree) {
    for (size_t k = 5; k < 30; ++k) {
        DBGHash graph(k);
        std::string bases("ACGT");
        std::string sequence;
        for (size_t i = 0; i < k * 3 / 2 + 1; ++i) {
            sequence.push_back(bases[(i * 199 % 163) % bases.size()]);
        }

        graph.add_sequence(sequence);

        for (size_t i = graph.first_edge(); i <= graph.last_edge(); ++i) {
            auto kmer = graph.get_node_kmer(i);
            auto label = graph.get_edge_label(i);
            if (kmer.front() != '$' && label != '$') {
                ASSERT_FALSE(graph.is_dummy_edge(kmer + label));
                ASSERT_TRUE(graph.has_the_only_outgoing_edge(i)) << kmer + label;
            }
        }
    }
}

TEST(DBGHash, SimplePathInDegree) {
    for (size_t k = 5; k < 30; ++k) {
        DBGHash graph(k);
        std::string bases("ACGT");
        std::string sequence;
        for (size_t i = 0; i < k * 3 / 2 + 1; ++i) {
            sequence.push_back(bases[(i * 199 % 163) % bases.size()]);
        }

        graph.add_sequence(sequence);

        for (size_t i = graph.first_edge(); i <= graph.last_edge(); ++i) {
            auto kmer = graph.get_node_kmer(i);
            auto label = graph.get_edge_label(i);
            if (kmer.front() != '$' && label != '$') {
                ASSERT_FALSE(graph.is_dummy_edge(kmer + label));
                ASSERT_TRUE(graph.has_the_only_outgoing_edge(i));
            }
        }
    }
}

TEST(DBGHash, TwoPathsDegree) {
    DBGHash graph(4);

    graph.add_sequence("AAAACGT");
    graph.add_sequence("AAAAGGT");

    ASSERT_FALSE(graph.has_the_only_outgoing_edge(graph.first_edge()))
        << graph.get_node_kmer(graph.first_edge())
        << graph.get_edge_label(graph.first_edge());
    ASSERT_TRUE(graph.has_the_only_incoming_edge(graph.first_edge()))
        << graph.get_node_kmer(graph.first_edge())
        << graph.get_edge_label(graph.first_edge());
    ASSERT_TRUE(graph.is_dummy_edge(
        graph.get_node_kmer(graph.first_edge())
            + graph.get_edge_label(graph.first_edge())
    )) << graph.get_node_kmer(graph.first_edge())
       << graph.get_edge_label(graph.first_edge());

    for (size_t i = graph.first_edge() + 5; i <= graph.last_edge(); ++i) {
        auto kmer = graph.get_node_kmer(i);
        auto label = graph.get_edge_label(i);
        if (kmer.front() != '$' && label != '$') {
            ASSERT_FALSE(graph.is_dummy_edge(kmer + label));
        }
        if (kmer[kmer.size() - 1] == 'A') {
            ASSERT_FALSE(graph.has_the_only_outgoing_edge(i))
                << i << ": " << kmer + label;
        } else {
            ASSERT_TRUE(graph.has_the_only_outgoing_edge(i))
                << i << ": " << kmer + label;
        }
        ASSERT_TRUE(graph.has_the_only_incoming_edge(i))
            << i << ": " << kmer + label;
    }
}

TEST(DBGHash, SimplePathForwardTraversal) {
    for (size_t k = 5; k < 30; ++k) {
        DBGHash graph(k);
        std::string bases("ACGT");
        std::string sequence;
        for (size_t i = 0; i < k * 3 / 2 + 1; ++i) {
            sequence.push_back(bases[(i * 199 % 163) % bases.size()]);
        }

        graph.add_sequence(sequence);

        for (size_t i = graph.first_edge(); i < graph.last_edge(); ++i) {
            auto label = graph.get_edge_label(i);
            ASSERT_EQ(i + 1, graph.next_edge(i, label))
                << graph.get_node_kmer(i) + graph.get_edge_label(i);
        }
    }
}

TEST(DBGHash, SimplePathBackwardTraversal) {
    for (size_t k = 5; k < 30; ++k) {
        DBGHash graph(k);
        std::string bases("ACGT");
        std::string sequence;
        for (size_t i = 0; i < k * 3 / 2 + 1; ++i) {
            sequence.push_back(bases[(i * 199 % 163) % bases.size()]);
        }

        graph.add_sequence(sequence);

        for (size_t i = graph.first_edge() + 1; i <= graph.last_edge(); ++i) {
            ASSERT_EQ(i - 1, graph.prev_edge(i));
        }
    }
}
