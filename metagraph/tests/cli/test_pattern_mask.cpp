#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/build.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "cli/pattern.hpp"
#include "cli/transform_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// The dummy-edge mask the pattern search requires (DESIGN-pattern-search.md §4), made the two
// ways a graph built without it gets it: `metagraph transform --mask-dummy`, which writes the
// .edgemask beside the graph, and --pattern-build-mask, which builds it in memory at load. Both
// must give exactly the mask `metagraph build --mask-dummy` writes: each test builds a graph with
// the real build command (masked), strips the mask from a copy, gives it back the other way and
// compares the bit vectors (and the files' bytes). Then the capabilities' `mask` and the
// mask_required message.

namespace {

using namespace mtg;
using namespace mtg::cli;
using mtg::graph::AnnotatedDBG;
using mtg::graph::CanonicalDBG;
using mtg::graph::DBGSuccinct;
using mtg::graph::DeBruijnGraph;
namespace fs = std::filesystem;

const size_t kK = 5;
// records whose dummy edges are not trivial: two starting alike with no k-mer before them
// (GAGAC, GAGTC: the source-dummy tree branches below $$GAG), one split by Ns (two islands,
// two starts and two ends), one shorter than k (no k-mer, no dummy), one of exactly k, one
// overlapping another's end, and starts that do have a k-mer before them (no dummy needed)
const std::vector<std::string> kRecords = {
    "ACGTACGGTTACCAGTAACGTTGCAAGT",
    "ACGTACGGTAACC",
    "TTGCANNNNACGGTCA",
    "ACG",
    "CAGTA",
    "GGTTACCAGTAACGTTGCAAGTCCCGGGAAATTT",
    "GAGACCTTAAGGC",
    "GAGTCATCCAAGG",
};

// a Config made from a command line, as main() makes it (the Config sets the swap path; the
// one the other tests use is restored)
std::unique_ptr<Config> make_config(std::vector<std::string> args) {
    const fs::path swap = utils::get_swap_path();
    args.insert(args.begin(), "metagraph");
    std::vector<char*> argv;
    for (std::string &a : args) {
        argv.push_back(a.data());
    }
    argv.push_back(nullptr);
    auto config = std::make_unique<Config>(static_cast<int>(args.size()), argv.data());
    utils::set_swap_path(swap);
    return config;
}

std::string read_file(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string rc(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

// A fresh directory with the records as FASTA
std::string make_dir(const std::string &name) {
    const std::string dir = test_dump_dir() + "/pattern_mask_" + name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    std::ofstream fasta(dir + "/records.fa");
    for (size_t i = 0; i < kRecords.size(); ++i) {
        fasta << ">r" << i << "\n" << kRecords[i] << "\n";
    }
    return dir;
}

// `metagraph build --mask-dummy --in-ram`: the reference mask
std::string build_masked(const std::string &dir, const std::string &mode,
                         const std::string &state) {
    const std::string base = dir + "/built";
    auto config = make_config({ "build", "--mask-dummy", "--in-ram", "-k", std::to_string(kK),
                                "--mode", mode, "--state", state, "-o", base,
                                dir + "/records.fa" });
    EXPECT_EQ(0, build_graph(config.get()));
    EXPECT_TRUE(fs::exists(base + ".edgemask"));
    return base + ".dbg";
}

// a copy of the built graph without its mask
std::string strip(const std::string &dir, const std::string &built) {
    const std::string stripped = dir + "/stripped.dbg";
    fs::copy_file(built, stripped, fs::copy_options::overwrite_existing);
    EXPECT_FALSE(fs::exists(dir + "/stripped.edgemask"));
    return stripped;
}

std::shared_ptr<DBGSuccinct> load(const std::string &path) {
    auto graph = std::make_shared<DBGSuccinct>(2);
    EXPECT_TRUE(graph->load(path)) << path;
    return graph;
}

void expect_same_mask(const DBGSuccinct &expected, const DBGSuccinct &actual) {
    ASSERT_NE(nullptr, expected.get_mask());
    ASSERT_NE(nullptr, actual.get_mask());
    ASSERT_EQ(expected.get_mask()->size(), actual.get_mask()->size());
    EXPECT_EQ(expected.get_mask()->num_set_bits(), actual.get_mask()->num_set_bits());
    EXPECT_TRUE(expected.get_mask()->to_vector() == actual.get_mask()->to_vector());
    EXPECT_EQ(expected.num_nodes(), actual.num_nodes());
}

const std::vector<std::string> kModes = { "basic", "canonical", "primary" };
const std::vector<std::string> kStates = { "stat", "fast", "small" };

TEST(PatternMask, TransformWritesTheMaskOfBuild) {
    for (const std::string &mode : kModes) {
        for (const std::string &state : kStates) {
            SCOPED_TRACE(mode + " " + state);
            const std::string dir = make_dir("transform_" + mode + "_" + state);
            const std::string built = build_masked(dir, mode, state);
            const std::string stripped = strip(dir, built);
            const std::string graph_bytes = read_file(stripped);

            ASSERT_EQ(0, transform_graph(make_config({ "transform", "--mask-dummy",
                                                       stripped }).get()));
            // only the mask was written: the graph file is as it was, nothing else appeared
            EXPECT_EQ(graph_bytes, read_file(stripped));
            std::set<std::string> files;
            for (const auto &entry : fs::directory_iterator(dir)) {
                files.insert(entry.path().filename().string());
            }
            EXPECT_EQ((std::set<std::string>{ "records.fa", "built.dbg", "built.edgemask",
                                              "stripped.dbg", "stripped.edgemask" }),
                      files);

            // the bits, and the very bytes, of the build's mask
            expect_same_mask(*load(built), *load(stripped));
            EXPECT_EQ(read_file(dir + "/built.edgemask"), read_file(dir + "/stripped.edgemask"));

            // an existing mask is not replaced unseen
            std::ofstream(dir + "/stripped.edgemask", std::ios::binary) << "not a mask";
            EXPECT_EQ(1, transform_graph(make_config({ "transform", "--mask-dummy",
                                                       stripped }).get()));
            EXPECT_EQ("not a mask", read_file(dir + "/stripped.edgemask"));
            // --force rebuilds it (without reading the broken one)
            ASSERT_EQ(0, transform_graph(make_config({ "transform", "--mask-dummy", "--force",
                                                       stripped }).get()));
            EXPECT_EQ(read_file(dir + "/built.edgemask"), read_file(dir + "/stripped.edgemask"));
            EXPECT_EQ(graph_bytes, read_file(stripped));
        }
    }
}

TEST(PatternMask, TransformWithMmapAndThreads) {
    const std::string dir = make_dir("transform_mmap");
    const std::string built = build_masked(dir, "basic", "stat");
    const std::string stripped = strip(dir, built);
    const bool mmap = utils::with_mmap();
    const size_t threads = get_num_threads();
    utils::set_mmap(true);
    set_num_threads(4);
    const int status = transform_graph(make_config({ "transform", "--mask-dummy",
                                                     stripped }).get());
    utils::set_mmap(mmap);
    set_num_threads(threads);
    ASSERT_EQ(0, status);
    EXPECT_EQ(read_file(dir + "/built.edgemask"), read_file(dir + "/stripped.edgemask"));
}

TEST(PatternMask, TransformRefusesWhatHasNoMask) {
    const std::string dir = make_dir("transform_refusals");
    // not a succinct graph
    std::ofstream(dir + "/graph.orhashdbg") << "x";
    EXPECT_EQ(1, transform_graph(make_config({ "transform", "--mask-dummy",
                                               dir + "/graph.orhashdbg" }).get()));
    // a missing graph: no mask is written for it
    EXPECT_EQ(1, transform_graph(make_config({ "transform", "--mask-dummy",
                                               dir + "/missing.dbg" }).get()));
    EXPECT_FALSE(fs::exists(dir + "/missing.edgemask"));
    EXPECT_FALSE(fs::exists(dir + "/graph.edgemask"));
}

TEST(PatternMask, BuildAtLoadEqualsTheMaskOfBuild) {
    for (const std::string &mode : kModes) {
        for (const std::string &state : kStates) {
            SCOPED_TRACE(mode + " " + state);
            const std::string dir = make_dir("load_" + mode + "_" + state);
            const std::string built = build_masked(dir, mode, state);
            const std::string stripped = strip(dir, built);

            std::shared_ptr<DeBruijnGraph> graph = load_critical_dbg(stripped);
            ASSERT_EQ(nullptr, dynamic_cast<DBGSuccinct&>(*graph).get_mask());
            EXPECT_FALSE(mask_built_at_load(*graph));
            build_mask_at_load(graph);
            EXPECT_TRUE(mask_built_at_load(*graph));
            expect_same_mask(*load(built), dynamic_cast<DBGSuccinct&>(*graph));
            // the served PRIMARY graph is the wrapper: it is recognised through it
            if (mode == "primary") {
                EXPECT_TRUE(mask_built_at_load(CanonicalDBG(graph)));
            }
            // nothing was written
            EXPECT_FALSE(fs::exists(dir + "/stripped.edgemask"));

            // a graph with its mask keeps it: nothing is built, and it is not "built at load"
            std::shared_ptr<DeBruijnGraph> with_file = load_critical_dbg(built);
            const auto *mask_before = dynamic_cast<DBGSuccinct&>(*with_file).get_mask();
            ASSERT_NE(nullptr, mask_before);
            build_mask_at_load(with_file);
            EXPECT_EQ(mask_before, dynamic_cast<DBGSuccinct&>(*with_file).get_mask());
            EXPECT_FALSE(mask_built_at_load(*with_file));
        }
    }
}

TEST(PatternMask, BuiltAtLoadIsNotInheritedByAnotherGraph) {
    const std::string dir = make_dir("load_lifetime");
    const std::string built = build_masked(dir, "basic", "stat");
    const std::string stripped = strip(dir, built);
    const DeBruijnGraph *address;
    {
        std::shared_ptr<DeBruijnGraph> graph = load_critical_dbg(stripped);
        build_mask_at_load(graph);
        address = graph.get();
        EXPECT_TRUE(mask_built_at_load(*graph));
    }
    // a graph loaded after the first was freed read its mask from the file. With the registry's
    // weak_ptr the expired entry pins the first graph's make_shared block until it is pruned,
    // so this graph is never at its address; at the same address only a registry of raw
    // pointers would be fooled, and whether a load lands there depends on the allocator (the
    // test below makes it land there)
    std::shared_ptr<DeBruijnGraph> graph = load_critical_dbg(built);
    EXPECT_FALSE(mask_built_at_load(*graph)) << (graph.get() == address ? "same address" : "");
}

// review of 2026-10-07, T3-03: the second graph is constructed at the first one's address by
// construction, whatever the allocator: one block of storage for both, their control blocks
// allocated apart (an expired registry entry pins the control block, never the storage). A
// registry that remembered addresses would take the second graph's mask for one built at load
TEST(PatternMask, BuiltAtLoadIsNotInheritedAtTheSameAddress) {
    const std::string dir = make_dir("load_address");
    const std::string built = build_masked(dir, "basic", "stat");
    const std::string stripped = strip(dir, built);
    alignas(DBGSuccinct) unsigned char storage[sizeof(DBGSuccinct)];
    auto destroy = [](DBGSuccinct *g) { g->~DBGSuccinct(); };
    {
        DBGSuccinct *first = new (storage) DBGSuccinct(2);
        std::shared_ptr<DeBruijnGraph> graph(first, destroy);
        ASSERT_TRUE(first->load(stripped));
        ASSERT_EQ(nullptr, first->get_mask());
        build_mask_at_load(graph);
        ASSERT_NE(nullptr, first->get_mask());
        EXPECT_TRUE(mask_built_at_load(*graph));
    }
    DBGSuccinct *second = new (storage) DBGSuccinct(2);
    std::shared_ptr<DeBruijnGraph> graph(second, destroy);
    ASSERT_TRUE(second->load(built));
    ASSERT_NE(nullptr, second->get_mask());
    ASSERT_EQ(static_cast<const void*>(storage),
              static_cast<const void*>(dynamic_cast<DBGSuccinct*>(graph.get())));
    EXPECT_FALSE(mask_built_at_load(*graph));
}

TEST(PatternMask, CountsOfTheMask) {
    for (const std::string &mode : kModes) {
        SCOPED_TRACE(mode);
        const std::string dir = make_dir("counts_" + mode);
        const std::string built = build_masked(dir, mode, "stat");
        auto graph = std::make_shared<DBGSuccinct>(2);
        ASSERT_TRUE(graph->load_without_mask(built));
        const auto &boss = graph->get_boss();
        const uint64_t source = boss.mark_source_dummy_edges(nullptr, 1);
        const uint64_t sink = boss.mark_sink_dummy_edges(nullptr);

        const DummyMaskCounts counts = mask_dummy_edges(graph.get(), 2);
        EXPECT_EQ(boss.num_edges(), counts.edges);
        // `metagraph stats --count-dummy`'s figures
        EXPECT_EQ(source, counts.source_dummy);
        EXPECT_EQ(sink, counts.sink_dummy);
        EXPECT_EQ(counts.edges - source - sink, counts.kmers);
        EXPECT_EQ(graph->num_nodes(), counts.kmers);
        EXPECT_GT(source, 1u);
        EXPECT_GT(sink, 0u);

        if (mode != "primary") {
            // the k-mers of the records' DNA4 islands (and their reverse complements)
            std::set<std::string> kmers;
            for (const std::string &record : kRecords) {
                for (size_t i = 0; i + kK <= record.size(); ++i) {
                    const std::string kmer = record.substr(i, kK);
                    if (kmer.find('N') != std::string::npos)
                        continue;
                    kmers.insert(kmer);
                    if (mode == "canonical")
                        kmers.insert(rc(kmer));
                }
            }
            EXPECT_EQ(kmers.size(), counts.kmers);
        }
    }
}

// The capabilities' mask and the refusal, through the loader the server and the CLI use
// (initialize_annotated_dbg with the Config of `metagraph pattern`)
TEST(PatternMask, CapabilitiesAndRefusal) {
    for (const std::string &mode : kModes) {
        SCOPED_TRACE(mode);
        const std::string dir = make_dir("caps_" + mode);
        const std::string built = build_masked(dir, mode, "stat");
        const std::string stripped = strip(dir, built);
        const std::string anno = dir + "/anno.column.annodbg";
        {
            // one label on every row: an annotation the loader accepts (no reads happen)
            auto graph = load(built);
            annot::ColumnCompressed<> annotation(graph->max_index());
            std::vector<uint64_t> rows(graph->max_index());
            std::iota(rows.begin(), rows.end(), 0);
            annotation.add_labels(rows, { "records" });
            annotation.serialize(anno);
        }
        std::ofstream(dir + "/request.json") << "{}";
        PatternLimits limits;
        limits.min_information_bits = 4;
        const std::string request = "{\"patterns\": [{\"dna\": \"ACGG\"}], \"mode\": \"count\","
                                    " \"scope\": \"any_offset\"}";
        auto served = [&](const std::string &graph, bool build_mask) {
            std::vector<std::string> args = { "pattern", "-i", graph, "-a", anno,
                                               dir + "/request.json" };
            if (build_mask)
                args.insert(args.begin() + 1, "--pattern-build-mask");
            return initialize_annotated_dbg(*make_config(args));
        };

        auto from_file = served(built, false);
        auto built_at_load = served(stripped, true);
        auto absent = served(stripped, false);
        // the flag on a graph with its file: the file's mask
        auto file_and_flag = served(built, true);

        auto mask_of = [&](const AnnotatedDBG &anno_graph) {
            return pattern_capabilities_json(&anno_graph, limits, false)["mask"].asString();
        };
        EXPECT_EQ("file", mask_of(*from_file));
        EXPECT_EQ("built_at_load", mask_of(*built_at_load));
        EXPECT_EQ("file", mask_of(*file_and_flag));
        EXPECT_EQ("absent", mask_of(*absent));
        EXPECT_EQ("mask_required", pattern_capabilities_json(absent.get(), limits, false)
                                           ["unavailable_reason"].asString());
        EXPECT_TRUE(pattern_capabilities_json(built_at_load.get(), limits, false)["available"]
                        .asBool());

        // the counts of a mask built at load are those of the file's
        const Json::Value json = parse_pattern_body(request);
        Json::Value a = process_pattern_request(json, *from_file, limits, "");
        Json::Value b = process_pattern_request(json, *built_at_load, limits, "");
        ASSERT_EQ(1u, a["patterns"].size());
        EXPECT_EQ(a["patterns"][0]["counts"], b["patterns"][0]["counts"]);
        EXPECT_EQ(a["patterns"][0]["work"], b["patterns"][0]["work"]);
        EXPECT_GT(a["patterns"][0]["counts"]["contexts"]["value"].asUInt64(), 0u);

        // refused without a mask, the message naming both remedies
        try {
            process_pattern_request(json, *absent, limits, "");
            ADD_FAILURE() << "a graph without its mask was searched";
        } catch (const PatternRefusal &e) {
            EXPECT_EQ(400, e.status());
            EXPECT_EQ("mask_required", e.code());
            const std::string message = e.what();
            EXPECT_NE(std::string::npos,
                      message.find("metagraph transform --mask-dummy <graph>.dbg")) << message;
            EXPECT_NE(std::string::npos, message.find("--pattern-build-mask")) << message;
            EXPECT_EQ(e.body()["error"].asString(), message);
        }
    }
}

} // namespace
