#include <fstream>
#include <map>
#include <set>

#include <gtest/gtest.h>

#include "../test_helpers.hpp"
#include "all/test_dbg_helpers.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include <include/builder/util.hpp>

namespace {

TEST(SSHashMinimizers, EmptyAndFinalTupleBuckets) {
    sshash::minimizers_tuples_iterator empty(nullptr, nullptr);
    EXPECT_FALSE(empty.has_next());
    const sshash::minimizer_tuple tuples[] = {
        { sshash::minimizer_info(1, 0, 0), 1 },
        { sshash::minimizer_info(1, 1, 0), 1 },
        { sshash::minimizer_info(2, 2, 0), 1 }
    };
    sshash::minimizers_tuples_iterator it(std::begin(tuples), std::end(tuples));
    ASSERT_TRUE(it.has_next());
    EXPECT_EQ(1u, it.minimizer());
    it.next();
    ASSERT_TRUE(it.has_next());
    EXPECT_EQ(2u, it.minimizer());
    it.next();
    EXPECT_FALSE(it.has_next());
}

#if ! _PROTEIN_GRAPH
TEST(SSHashPath, ReverseComplementPreservesNodeSpellings) {
    using namespace mtg::graph;
    for (auto mode : { DBGSSHash::BASIC, DBGSSHash::CANONICAL }) {
        SCOPED_TRACE(mode);
        auto graph = mtg::test::build_graph<DBGSSHash>(3, { "AACGTA" }, mode);
        for (std::string sequence : { "AACGTA", "AACNTAC" }) {
            auto path = map_to_nodes_sequentially(*graph, sequence);
            auto expected_sequence = sequence;
            reverse_complement(expected_sequence.begin(), expected_sequence.end());
            auto expected_path = map_to_nodes_sequentially(*graph, expected_sequence);
            reverse_complement_seq_path(*graph, sequence, path);
            EXPECT_EQ(expected_sequence, sequence);
            EXPECT_EQ(expected_path, path);
        }
    }
}

TEST(SSHashMinimizers, EnumeratesAllBucketSizes) {
    using mtg::graph::DBGSSHash;
    constexpr size_t k = 7;
    for (auto mode : { DBGSSHash::BASIC, DBGSSHash::CANONICAL }) {
        SCOPED_TRACE(mode);
        const std::string filename = test_dump_dir() + "/minimizer_positions.fa";
        {
            std::ofstream out(filename);
            for (size_t value = 0; value < (1u << (2 * k)); ++value) {
                std::string sequence(k, 'A');
                for (size_t i = 0; i < k; ++i)
                    sequence[i] = "ACGT"[(value >> (2 * i)) & 3];
                auto rc = sequence;
                reverse_complement(rc.begin(), rc.end());
                if (mode == DBGSSHash::BASIC || sequence <= rc)
                    out << ">\n" << sequence << '\n';
            }
        }
        // Select m=3 while retaining enough distinct k-mers to exercise
        // singleton, ordinary and skew-index buckets.
        DBGSSHash graph(filename, k, mode, 8, test_dump_dir());
        std::visit([&](const auto &dict) {
            using Kmer = typename std::decay_t<decltype(dict)>::kmer_type;
            ASSERT_EQ(3u, dict.m());
            std::map<uint64_t, std::set<uint64_t>> expected;
            for (uint64_t id = 0; id < dict.num_kmers(); ++id) {
                std::string sequence(k, 'A');
                dict.access(id, sequence.data());
                auto word = sshash::util::string_to_uint_kmer<Kmer>(sequence.data(), k);
                auto mini = sshash::util::compute_minimizer(word, k, dict.m(), dict.hasher());
                if (mode == DBGSSHash::CANONICAL) {
                    word.reverse_complement_inplace(k);
                    auto rc = sshash::util::compute_minimizer(word, k, dict.m(), dict.hasher());
                    if (rc.minimizer < mini.minimizer) {
                        mini = rc;
                        mini.pos_in_kmer = k - dict.m() - mini.pos_in_kmer;
                    }
                }
                auto [string_id, offset] = dict.strings_offsets().id_to_offset(id, k);
                std::ignore = string_id;
                expected[mini.minimizer].insert(offset + mini.pos_in_kmer);
            }
            bool singleton = false, ordinary = false, heavy = false;
            for (uint64_t mini = 0; mini < (1u << (2 * dict.m())); ++mini) {
                const auto &wanted = expected[mini];
                singleton |= wanted.size() == 1;
                ordinary |= wanted.size() > 1 && wanted.size() <= (1u << sshash::constants::min_l);
                heavy |= wanted.size() > (1u << sshash::constants::min_l);
                std::set<uint64_t> actual;
                dict.call_minimizer_positions(mini, [&](uint64_t offset) {
                    EXPECT_TRUE(actual.insert(offset).second);
                });
                EXPECT_EQ(wanted, actual) << "minimizer=" << mini;
            }
            EXPECT_TRUE(singleton);
            EXPECT_TRUE(ordinary);
            EXPECT_TRUE(heavy);
        }, graph.data());
    }
}
#endif

} // namespace
