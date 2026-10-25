#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <random>


#include <gtest/gtest.h>

#include "../test_helpers.hpp"
#include "all/test_dbg_helpers.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include <include/builder/util.hpp>
#include <kmer_iterator.hpp>
#include <streaming_query.hpp>
// Lookup is a template defined outside dictionary.hpp. Include its definition
// here rather than relying on instantiations emitted by another object file.
#include <src/dictionary.cpp>
#include <sdsl/memory_management.hpp>

namespace {

template <typename Word>
void check_reverse_refill() {
    using Kmer = mtg::graph::DBGSSHash::kmer_t<Word>;
    constexpr size_t k = Kmer::max_k;
    std::mt19937 random(672);
    std::string sequence(400, 'A');
    bits::bit_vector::builder builder;
    for (char &c : sequence) {
        c = Kmer::alphabet[random() % Kmer::alphabet_size];
        builder.append_bits(Kmer::char_to_uint(c), Kmer::bits_per_char);
    }
    // SSHash normally supplies padded storage for word-sized reads.
    for (size_t i = 0; i < Kmer::uint_kmer_bits / 64 + 1; ++i)
        builder.append_bits(0, 64);
    bits::bit_vector strings;
    builder.build(strings);
    sshash::kmer_iterator<Kmer, bits::bit_vector> it(strings, k,
                                                   sequence.size() * Kmer::bits_per_char);
    for (size_t end = sequence.size(); end >= k; --end) {
        auto expected = sshash::util::string_to_uint_kmer<Kmer>(sequence.data() + end - k, k);
        ASSERT_TRUE(expected == it.get_reverse()) << "end=" << end << " k=" << k;
        it.next_reverse();
    }
}

TEST(SSHashIterator, ReverseRefill64) { check_reverse_refill<uint64_t>(); }
TEST(SSHashIterator, ReverseRefill128) { check_reverse_refill<__uint128_t>(); }
TEST(SSHashIterator, ReverseRefill256) { check_reverse_refill<sdsl::uint256_t>(); }

#if _DNA5_GRAPH
TEST(SSHashFixtures, HighKPreservesStrandsAndBranches) {
    using namespace mtg::graph;
    constexpr size_t k = 87; // Beyond DBGHashFast's DNA5 packing limit.
    std::mt19937 random(672);
    std::string sequence(200, 'A');
    for (char &c : sequence)
        c = std::string_view("ACGT")[random() % 4];
    auto branch = sequence;
    branch[110] = branch[110] == 'A' ? 'C' : 'A';
    auto overlapping_reverse = sequence.substr(20, 140);
    reverse_complement(overlapping_reverse.begin(), overlapping_reverse.end());
    const std::vector<std::string> sequences { sequence, branch, overlapping_reverse };
    for (auto mode : { DBGSSHash::BASIC, DBGSSHash::CANONICAL, DBGSSHash::PRIMARY }) {
        SCOPED_TRACE(mode);
        std::set<std::string> expected;
        for (const auto &input : sequences) {
            for (size_t i = 0; i + k <= input.size(); ++i) {
                auto kmer = input.substr(i, k);
                expected.insert(kmer);
                if (mode != DBGSSHash::BASIC) {
                    reverse_complement(kmer.begin(), kmer.end());
                    expected.insert(kmer);
                }
            }
        }
        auto graph = mtg::test::build_graph<DBGSSHash>(k, sequences, mode);
        std::set<std::string> actual;
        graph->call_nodes([&](auto node) { actual.insert(graph->get_node_sequence(node)); });
        EXPECT_EQ(expected, actual);
        EXPECT_EQ(expected.size(), graph->num_nodes());
        for (const auto &kmer : expected) {
            const auto node = graph->kmer_to_node(kmer);
            ASSERT_NE(DeBruijnGraph::npos, node);
            for (char c : std::string_view("ACGT")) {
                const auto next = kmer.substr(1) + c;
                EXPECT_EQ(expected.contains(next), graph->traverse(node, c) != DeBruijnGraph::npos);
            }
        }
    }
}
#endif

TEST(SSHashFormat, AcceptsSameMajorVersion) {
    using mtg::graph::DBGSSHash;
    auto graph = mtg::test::build_graph<DBGSSHash>(3, { "ACGT" }, DBGSSHash::BASIC);
    std::stringstream bytes;
    static_cast<const DBGSSHash&>(*graph).serialize(bytes);
    bytes.seekp(2 * sizeof(size_t) + sizeof(DBGSSHash::Mode) + 1);
    bytes.put(sshash::constants::current_version_number::y + 1);
    DBGSSHash loaded(3);
    ASSERT_TRUE(loaded.load(bytes));
    ASSERT_EQ(graph->num_nodes(), loaded.num_nodes());
    for (size_t node = 1; node <= graph->num_nodes(); ++node)
        EXPECT_EQ(graph->get_node_sequence(node), loaded.get_node_sequence(node));
}

TEST(SSHashFormat, RejectsDifferentMajorVersion) {
    using mtg::graph::DBGSSHash;
    auto graph = mtg::test::build_graph<DBGSSHash>(3, { "ACGT" }, DBGSSHash::BASIC);
    std::stringstream bytes;
    static_cast<const DBGSSHash&>(*graph).serialize(bytes);
    // The wrapper's node count, k, and mode precede the dictionary version.
    bytes.seekp(2 * sizeof(size_t) + sizeof(DBGSSHash::Mode));
    bytes.put(sshash::constants::current_version_number::x + 1);
    DBGSSHash loaded(3);
    EXPECT_FALSE(loaded.load(bytes));
    EXPECT_EQ(0u, loaded.num_nodes());
}

TEST(SSHashFormat, RejectsV5BeforeReadingPayload) {
    using mtg::graph::DBGSSHash;
    auto graph = mtg::test::build_graph<DBGSSHash>(3, { "ACGT" }, DBGSSHash::BASIC);
    std::stringstream bytes;
    static_cast<const DBGSSHash&>(*graph).serialize(bytes);
    bytes.seekp(2 * sizeof(size_t) + sizeof(DBGSSHash::Mode));
    bytes.put(5);
    // Only the header and version are needed to reject an obsolete index.
    std::string header = bytes.str().substr(0, 2 * sizeof(size_t) + sizeof(DBGSSHash::Mode) + 3);
    const auto filename = test_dump_dir() + "/sshash_v5_header.sshashdbg";
    { std::ofstream out(filename, std::ios::binary); out.write(header.data(), header.size()); }
    DBGSSHash loaded(3);
    std::ifstream file(filename, std::ios::binary);
    EXPECT_FALSE(loaded.load(file));
    sdsl::mmap_ifstream mapped(filename, std::ios::binary);
    EXPECT_FALSE(loaded.load(mapped));
}

TEST(SSHashFormat, FailedReloadPreservesGraph) {
    using mtg::graph::DBGSSHash;
    auto original = mtg::test::build_graph<DBGSSHash>(3, { "ACGT" }, DBGSSHash::BASIC);
    std::stringstream bytes;
    static_cast<const DBGSSHash&>(*original).serialize(bytes);
    const auto serialized = bytes.str();
    DBGSSHash loaded(3);
    ASSERT_TRUE(loaded.load(bytes));

    const size_t header_size = 2 * sizeof(size_t) + sizeof(DBGSSHash::Mode);
    auto old_version = serialized;
    old_version[header_size] = 5;
    for (const auto &bad : { std::string(), old_version,
                            serialized.substr(0, header_size + 3 + sizeof(uint64_t)),
                            serialized.substr(0, serialized.size() - 1) }) {
        std::stringstream input(bad);
        EXPECT_FALSE(loaded.load(input));
        std::stringstream after;
        loaded.serialize(after);
        EXPECT_EQ(serialized, after.str());
        EXPECT_TRUE(loaded.find("ACGT"));
    }
}

TEST(SSHashFormat, LoadPreservesStreamExceptionMask) {
    using mtg::graph::DBGSSHash;
    auto graph = mtg::test::build_graph<DBGSSHash>(3, { "ACGT" }, DBGSSHash::BASIC);
    std::stringstream bytes;
    static_cast<const DBGSSHash&>(*graph).serialize(bytes);
    for (auto mask : { std::ios::goodbit, std::ios::badbit, std::ios::failbit,
                       std::ios::eofbit | std::ios::failbit | std::ios::badbit }) {
        for (bool truncated : { false, true }) {
            SCOPED_TRACE(mask);
            SCOPED_TRACE(truncated);
            std::stringstream input(truncated ? bytes.str().substr(0, 1) : bytes.str());
            input.exceptions(mask);
            DBGSSHash loaded(3);
            EXPECT_EQ(!truncated, loaded.load(input));
            EXPECT_EQ(mask, input.exceptions());
            EXPECT_EQ(truncated, input.fail());
        }
    }
}

TEST(SSHashLookup, ForwardOrientationsAndRoundTrip) {
    using mtg::graph::DBGSSHash;
#if _PROTEIN_GRAPH
    const std::string sequence = "MKWVTFISRLDPQEAHNYCG";
#else
    const std::string sequence = "ACGT";
#endif
    auto graph = mtg::test::build_graph<DBGSSHash>(3, { sequence }, DBGSSHash::BASIC);
    const auto filename = test_dump_dir() + "/sshash_v6_roundtrip.sshashdbg";
    graph->serialize(filename);
    auto check = [&](const DBGSSHash &candidate) {
        std::visit([&](const auto &dict) {
            for (size_t i = 0; i + 3 <= sequence.size(); ++i) {
                auto hit = dict.lookup(sequence.data() + i, false);
                ASSERT_NE(sshash::constants::invalid_uint64, hit.kmer_id);
                EXPECT_EQ(sshash::constants::forward_orientation, hit.kmer_orientation);
                EXPECT_EQ(sequence.substr(i, 3),
                          candidate.get_node_sequence(DBGSSHash::sshash_to_graph_index(hit.kmer_id)));
            }
#if _PROTEIN_GRAPH
            sshash::streaming_query<std::decay_t<decltype(dict)>> stream(&dict);
            for (size_t i = 0; i + 3 <= sequence.size(); ++i)
                EXPECT_TRUE(sshash::equal_lookup_result(dict.lookup(sequence.data() + i, false),
                                                        stream.lookup(sequence.data() + i)));
#endif
        }, candidate.data());
    };
    check(static_cast<const DBGSSHash&>(*graph));
    DBGSSHash loaded(3);
    std::ifstream file(filename, std::ios::binary);
    ASSERT_TRUE(loaded.load(file));
    check(loaded);
    sdsl::mmap_ifstream mapped(filename, std::ios::binary);
    ASSERT_TRUE(loaded.load(mapped));
    check(loaded);
}

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

TEST(SSHashSeeds, AmbiguousShortQueries) {
    using namespace mtg::graph;
    using namespace mtg::graph::align_redone;
    const std::string reference =
        "AACCGACAGAACTGATGGCAGGGCCTAGTCTCCGAACTTCTGTTAAGACAGGCGAAAGTTAACC";
    const std::string filename = test_dump_dir() + "/ambiguous_short_seeds.fa";
    { std::ofstream out(filename); out << ">reference\n" << reference << '\n'; }
#if _PROTEIN_GRAPH
    const char invalid = '*';
    const std::vector<DBGSSHash::Mode> modes { DBGSSHash::BASIC };
#else
    const char invalid = 'N';
    const std::vector<DBGSSHash::Mode> modes { DBGSSHash::BASIC, DBGSSHash::CANONICAL };
#endif
    for (auto mode : modes) {
        SCOPED_TRACE(mode);
        // Choose m=3 so an invalid character can lie inside a seed window.
        DBGSSHash graph(filename, 7, mode, 8, test_dump_dir());
        std::visit([](const auto &dict) { ASSERT_EQ(3u, dict.m()); }, graph.data());
        DBGAlignerConfig config;
        config.min_seed_length = 3;
        config.max_seed_length = graph.get_k();
        size_t valid_flank_seeds = 0;
        for (size_t length = 3; length < graph.get_k(); ++length) {
            for (size_t start = 0; start + length <= reference.size(); ++start) {
                for (size_t ambiguous = 0; ambiguous < length; ++ambiguous) {
                    std::string sequence = reference.substr(start, length);
                    sequence[ambiguous] = invalid;
                    SCOPED_TRACE(sequence);
                    Query query(graph, sequence);
                    for (const auto &anchor : ExactSeeder(query, config).get_anchors()) {
                        ASSERT_EQ(1u, anchor.get_path().size());
                        EXPECT_EQ(std::string::npos, anchor.get_seed().find(invalid));
                        EXPECT_EQ(anchor.get_seed(),
                                  graph.get_node_sequence(anchor.get_path().front())
                                      .substr(0, anchor.get_seed().size()));
                        ++valid_flank_seeds;
                    }
                }
            }
        }
        EXPECT_GT(valid_flank_seeds, 0u);
        Query ambiguous(graph, std::string(6, invalid));
        EXPECT_TRUE(ExactSeeder(ambiguous, config).get_anchors().empty());
        DBGSSHash empty(7, mode);
        Query empty_query(empty, "ACGTT");
        EXPECT_TRUE(ExactSeeder(empty_query, config).get_anchors().empty());
    }
}

TEST(SSHashPath, FinalNodeKeepsItsUnitigId) {
    using namespace mtg::graph;
    const std::string filename = test_dump_dir() + "/final_unitig_node.fa";
    { std::ofstream out(filename); out << ">reference\nACGTA\n"; }
    DBGSSHash graph(filename, 3, DBGSSHash::BASIC, 5, test_dump_dir());
    const auto path = map_to_nodes_sequentially(graph, "ACGTA");
    ASSERT_EQ(3u, path.size());
    ASSERT_EQ(graph.num_nodes(), path.back());
    for (auto node : path)
        EXPECT_EQ(graph.get_sequence_ids(path.front()), graph.get_sequence_ids(node));
}

TEST(SSHashPath, LoadRepairsLegacyFinalPredecessor) {
    using namespace mtg::graph;
    const std::string filename = test_dump_dir() + "/legacy_final_unitig_node.fa";
    { std::ofstream out(filename); out << ">reference\nACGTA\n"; }
    DBGSSHash graph(filename, 3, DBGSSHash::BASIC, 5, test_dump_dir());
    ASSERT_EQ(3u, graph.num_nodes());
    ASSERT_EQ("ACG", graph.get_node_sequence(1));
    ASSERT_EQ("CGT", graph.get_node_sequence(2));
    ASSERT_EQ("GTA", graph.get_node_sequence(3));

    std::stringstream bytes;
    essentials::generic_saver saver(bytes);
    saver.visit(graph.num_nodes());
    saver.visit(graph.get_k());
    saver.visit(graph.get_mode());
    std::visit([&](const auto &dict) { saver.visit(dict); }, graph.data());
    const auto cache_offset = bytes.tellp();
    bit_vector_smart({ false, true, true, false }).serialize(bytes);
    // Older graphs left the last predecessor bit unset even on a linear path.
    bit_vector_smart({ false, false, true, false }).serialize(bytes);

    DBGSSHash loaded(3);
    ASSERT_TRUE(loaded.load(bytes));
    for (auto node : map_to_nodes_sequentially(loaded, "ACGTA"))
        EXPECT_EQ(loaded.get_sequence_ids(1), loaded.get_sequence_ids(node));
    EXPECT_EQ(2u, loaded.traverse_back(3, 'C'));
    std::stringstream repaired;
    loaded.serialize(repaired);
    repaired.seekg(cache_offset);
    bit_vector_smart successors, predecessors;
    ASSERT_TRUE(successors.load(repaired));
    ASSERT_TRUE(predecessors.load(repaired));
    ASSERT_EQ(4u, successors.size());
    ASSERT_EQ(4u, predecessors.size());
    // Check persisted bits, excluding unused padding in the last storage word.
    EXPECT_EQ(0b0110u, successors.get_int(0, 4));
    EXPECT_EQ(0b1100u, predecessors.get_int(0, 4));
    repaired.seekg(0);
    DBGSSHash roundtrip(3);
    ASSERT_TRUE(roundtrip.load(repaired));
    for (auto node : map_to_nodes_sequentially(roundtrip, "ACGTA")) {
        EXPECT_EQ(graph.get_node_sequence(node), roundtrip.get_node_sequence(node));
        EXPECT_EQ(graph.get_sequence_ids(node), roundtrip.get_sequence_ids(node));
    }
}

#if ! _PROTEIN_GRAPH
TEST(SSHashPath, FinalNodeCanLeadToTheReverseStrand) {
    using namespace mtg::graph;
    const std::string filename = test_dump_dir() + "/final_unitig_reverse.fa";
    { std::ofstream out(filename); out << ">reference\nAAT\n"; }
    DBGSSHash graph(filename, 3, DBGSSHash::CANONICAL, 3, test_dump_dir());
    const auto forward = graph.kmer_to_node("AAT");
    const auto reverse = graph.kmer_to_node("ATT");
    ASSERT_EQ(1u, forward);
    ASSERT_EQ(2u, reverse);
    EXPECT_EQ(reverse, graph.traverse(forward, 'T'));
    std::vector<DeBruijnGraph::node_index> outgoing;
    graph.call_outgoing_kmers(forward, [&](auto node, char c) {
        EXPECT_EQ('T', c);
        outgoing.push_back(node);
    });
    EXPECT_EQ(std::vector<DeBruijnGraph::node_index>{ reverse }, outgoing);
}

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
            // Canonical minimizers put reverse-complement pairs in the same
            // bucket. Retain one kmer of the least frequent bucket to cover a
            // singleton even in BASIC mode; keep all other pairs for skew tests.
            using Word = DBGSSHash::kmer_t<uint64_t>;
            std::map<uint64_t, std::vector<std::string>> buckets;
            for (size_t value = 0; value < (1u << (2 * k)); ++value) {
                std::string sequence(k, 'A');
                for (size_t i = 0; i < k; ++i)
                    sequence[i] = "ACGT"[(value >> (2 * i)) & 3];
                auto rc = sequence;
                reverse_complement(rc.begin(), rc.end());
                if (mode == DBGSSHash::BASIC || sequence <= rc) {
                    auto word = sshash::util::string_to_uint_kmer<Word>(sequence.data(), k);
                    auto mini = sshash::util::compute_minimizer(word, k, 3,
                                             sshash::hasher_type(sshash::constants::seed));
                    buckets[mini.minimizer].push_back(sequence);
                }
            }
            auto smallest = std::min_element(buckets.begin(), buckets.end(),
                [](const auto &a, const auto &b) { return a.second.size() < b.second.size(); });
            smallest->second.resize(1);
            for (const auto &[mini, sequences] : buckets)
                for (const auto &sequence : sequences)
                    out << ">\n" << sequence << '\n';
        }
        // Select m=3 while retaining enough distinct k-mers to exercise
        // singleton, ordinary and skew-index buckets.
        DBGSSHash graph(filename, k, mode, 8, test_dump_dir());
        std::visit([&](const auto &dict) {
            using Kmer = typename std::decay_t<decltype(dict)>::kmer_type;
            ASSERT_EQ(3u, dict.m());
            std::map<uint64_t, std::set<uint64_t>> expected;
            for (uint64_t id = 0; id < dict.num_kmers(); ++id) {
                auto sequence = graph.get_node_sequence(DBGSSHash::sshash_to_graph_index(id));
                // Includes both orientations of every DNA k-mer in BASIC mode,
                // exercising duplicate canonical keys in the heavy skew index.
                EXPECT_EQ(id, dict.lookup(sequence.data(), false).kmer_id) << sequence;
                auto word = sshash::util::string_to_uint_kmer<Kmer>(sequence.data(), k);
                auto mini = sshash::util::compute_minimizer(word, k, dict.m(), dict.hasher());
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
