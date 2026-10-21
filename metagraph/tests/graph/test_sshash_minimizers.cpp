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
