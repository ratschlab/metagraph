#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "../graph/all/test_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/augment.hpp"
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


// The dummy-edge mask that makes the pattern search's counts exact (DESIGN-pattern-search.md
// §4), made the two ways a graph built without it gets it: `metagraph transform
// --mask-dummy`, which writes the .edgemask beside the graph, and --pattern-build-mask, which
// builds it in memory at load. Both must give exactly the mask `metagraph build --mask-dummy`
// writes: each test builds a graph with the real build command (masked), strips the mask from
// a copy, gives it back the other way and compares the bit vectors (and the files' bytes).
// Then the capabilities' `mask` and `counting`, and a graph without its mask, which is served
// (its counts upper bounds with an estimate from its sampled dummy fraction).

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

// The second graph is constructed at the first one's address by construction, whatever the
// allocator: one block of storage for both, their control blocks allocated apart (an expired
// registry entry pins the control block, never the storage). A registry that remembered
// addresses would take the second graph's mask for one built at load
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

// The capabilities' mask and counting, through the loader the server and the CLI use
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
        // a graph without its mask is served, its counts upper bounds
        const Json::Value absent_caps = pattern_capabilities_json(absent.get(), limits, false);
        EXPECT_TRUE(absent_caps["available"].asBool());
        EXPECT_TRUE(absent_caps["unavailable_reason"].isNull());
        EXPECT_EQ("upper_bound", absent_caps["counting"].asString());
        EXPECT_TRUE(pattern_capabilities_json(built_at_load.get(), limits, false)["available"]
                        .asBool());
        for (const auto *masked : { from_file.get(), built_at_load.get(), file_and_flag.get() }) {
            const Json::Value caps = pattern_capabilities_json(masked, limits, false);
            EXPECT_EQ("exact", caps["counting"].asString());
            EXPECT_TRUE(caps["dummy_fraction"].isNull());
        }
        // the loader sampled its dummy fraction in its thread: the one the engine samples
        const auto &base = dynamic_cast<const DBGSuccinct&>(
                mode == "primary"
                    ? dynamic_cast<const CanonicalDBG&>(absent->get_graph()).get_graph()
                    : absent->get_graph());
        const auto loaded = dummy_fraction(*absent);
        ASSERT_TRUE(loaded);
        const auto sampled = mtg::graph::pattern::sample_real_fraction(base);
        EXPECT_EQ(sampled.real, loaded->real);
        EXPECT_EQ(sampled.samples, loaded->samples);
        EXPECT_EQ(dummy_fraction_json(sampled), absent_caps["dummy_fraction"]);
        EXPECT_FALSE(dummy_fraction(*from_file));
        EXPECT_FALSE(dummy_fraction(*built_at_load));

        // the counts of a mask built at load are those of the file's
        const Json::Value json = parse_pattern_body(request);
        Json::Value a = process_pattern_request(json, *from_file, limits, "");
        Json::Value b = process_pattern_request(json, *built_at_load, limits, "");
        ASSERT_EQ(1u, a["patterns"].size());
        EXPECT_EQ(a["patterns"][0]["counts"], b["patterns"][0]["counts"]);
        EXPECT_EQ(a["patterns"][0]["work"], b["patterns"][0]["work"]);
        EXPECT_GT(a["patterns"][0]["counts"]["contexts"]["value"].asUInt64(), 0u);

        // answered without a mask: an upper bound holding the exact count, its estimate, and
        // the answer saying so; the masked answers say nothing of it
        const Json::Value c = process_pattern_request(json, *absent, limits, "");
        const Json::Value &exact = a["patterns"][0]["counts"]["contexts"];
        const Json::Value &bound = c["patterns"][0]["counts"]["contexts"];
        EXPECT_EQ("upper_bound", c["index"]["counting"].asString());
        EXPECT_EQ(absent_caps["dummy_fraction"], c["index"]["dummy_fraction"]);
        EXPECT_FALSE(a["index"].isMember("counting"));
        if (bound["relation"].asString() == "exact") {
            EXPECT_EQ(exact["value"], bound["value"]);
        } else {
            ASSERT_EQ("bounds", bound["relation"].asString()) << bound;
            EXPECT_LE(bound["lower"].asUInt64(), exact["value"].asUInt64());
            EXPECT_GE(bound["upper"].asUInt64(), exact["value"].asUInt64());
            EXPECT_TRUE(bound.isMember("estimate"));
        }
    }
}

// A mask whose W = $ edge is valid is refused (SPEC-pattern-search §10.2): such a mask, as
// DBGSuccinct::add_sequence writes on a masked graph, would make the pattern search claim a
// too-large count exact (CG: exact 3 where the graph holds 2). It is found once at load (O(W =
// $ edges)) and refused, mask_invalid, in the capabilities and with a 400 naming the remedy;
// `metagraph extend` rebuilds the mask of its output, and removes a mask left beside it that is
// not its graph's
TEST(PatternMask, MaskWithAValidSentinelIsRefused) {
    const std::string dir = test_dump_dir() + "/pattern_mask_invalid";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::vector<std::string> records = { "ACGTTGCA", "ACGTTGAC" };
    std::ofstream(dir + "/s1.fa") << ">r1\n" << records[0] << "\n";
    std::ofstream(dir + "/s2.fa") << ">r2\n" << records[1] << "\n";
    const std::string base = dir + "/g1";
    ASSERT_EQ(0, build_graph(make_config({ "build", "--mask-dummy", "--in-ram", "-k", "4",
                                           "-o", base, dir + "/s1.fa" }).get()));
    EXPECT_EQ(0u, load(base + ".dbg")->count_valid_sentinel_edges());

    // the mask add_sequence writes: every inserted edge valid, the new sink dummy GAC$ too
    {
        auto extended = load(base + ".dbg");
        extended->switch_state(mtg::graph::boss::BOSS::State::DYN);
        extended->add_sequence(records[1]);
        ASSERT_GT(extended->count_valid_sentinel_edges(), 0u);
        extended->serialize(dir + "/bad");
        ASSERT_TRUE(fs::exists(dir + "/bad.edgemask"));
        // the same graph re-masked: none
        extended->mask_dummy_kmers(1, false);
        EXPECT_EQ(0u, extended->count_valid_sentinel_edges());
    }

    // `metagraph extend` writes the mask transform --mask-dummy would build for its output
    ASSERT_EQ(0, augment_graph(make_config({ "extend", "-i", base + ".dbg", "-o", dir + "/ext",
                                             dir + "/s2.fa" }).get()));
    {
        auto ext = load(dir + "/ext.dbg");
        ASSERT_NE(nullptr, ext->get_mask());
        EXPECT_EQ(0u, ext->count_valid_sentinel_edges());
        auto remasked = std::make_shared<DBGSuccinct>(2);
        ASSERT_TRUE(remasked->load_without_mask(dir + "/ext.dbg"));
        mask_dummy_edges(remasked.get(), 1);
        expect_same_mask(*remasked, *ext);
        // the k-mers of the records, no dummy
        std::set<std::string> kmers;
        for (const std::string &r : records) {
            for (size_t i = 0; i + 4 <= r.size(); ++i) {
                kmers.insert(r.substr(i, 4));
            }
        }
        EXPECT_EQ(kmers.size(), ext->num_nodes());
    }
    // an input without a mask: a mask left beside the output by another graph is removed
    fs::copy_file(base + ".dbg", dir + "/nomask.dbg");
    fs::copy_file(base + ".edgemask", dir + "/ext2.edgemask");
    ASSERT_EQ(0, augment_graph(make_config({ "extend", "-i", dir + "/nomask.dbg", "-o",
                                             dir + "/ext2", dir + "/s2.fa" }).get()));
    EXPECT_TRUE(fs::exists(dir + "/ext2.dbg"));
    EXPECT_FALSE(fs::exists(dir + "/ext2.edgemask"));

    // served through the loader of the server and the CLI
    auto annotate_all = [&](const std::string &graph_path, const std::string &anno) {
        auto graph = load(graph_path);
        annot::ColumnCompressed<> annotation(graph->max_index());
        std::vector<uint64_t> rows(graph->max_index());
        std::iota(rows.begin(), rows.end(), 0);
        annotation.add_labels(rows, { "records" });
        annotation.serialize(anno);
    };
    annotate_all(dir + "/bad.dbg", dir + "/bad.column.annodbg");
    annotate_all(dir + "/ext.dbg", dir + "/ext.column.annodbg");
    std::ofstream(dir + "/request.json") << "{}";
    auto served = [&](const std::string &name) {
        return initialize_annotated_dbg(*make_config({ "pattern", "-i", dir + "/" + name + ".dbg",
                                                       "-a", dir + "/" + name + ".column.annodbg",
                                                       dir + "/request.json" }));
    };
    auto bad = served("bad");
    auto good = served("ext");
    EXPECT_TRUE(mask_invalid_at_load(bad->get_graph()));
    EXPECT_FALSE(mask_invalid_at_load(good->get_graph()));

    PatternLimits limits;
    limits.min_information_bits = 0;
    const Json::Value bad_caps = pattern_capabilities_json(bad.get(), limits, false);
    EXPECT_FALSE(bad_caps["available"].asBool());
    EXPECT_EQ("mask_invalid", bad_caps["unavailable_reason"].asString());
    EXPECT_EQ("file", bad_caps["mask"].asString());
    EXPECT_EQ("basic", bad_caps["graph_mode"].asString());
    EXPECT_TRUE(pattern_capabilities_json(good.get(), limits, false)["available"].asBool());

    const Json::Value json = parse_pattern_body(
            "{\"patterns\": [{\"dna\": \"CG\"}], \"mode\": \"count\", \"scope\": \"any_offset\"}");
    try {
        process_pattern_request(json, *bad, limits, "");
        ADD_FAILURE() << "a graph with an invalid mask was searched";
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(400, e.status());
        EXPECT_EQ("mask_invalid", e.code());
        const std::string message = e.what();
        EXPECT_NE(std::string::npos,
                  message.find("metagraph transform --mask-dummy --force <graph>.dbg"))
                << message;
    }
    // the extended graph's own mask answers the truth: CG (its own reverse complement) in
    // ACGT and CGTT
    const Json::Value answer = process_pattern_request(json, *good, limits, "");
    ASSERT_EQ(1u, answer["patterns"].size());
    EXPECT_EQ("exact", answer["patterns"][0]["counts"]["contexts"]["relation"].asString());
    EXPECT_EQ(2u, answer["patterns"][0]["counts"]["contexts"]["value"].asUInt64());
}


// ---------------------------------------------------------------- without the mask (#16)

// The Wilson interval at its extremes: every sample real (the interval's top is 1), and its
// JSON {value, interval, samples, source: "sampled"}
TEST(PatternMaskUnmasked, DummyFractionJson) {
    DummyFraction f;
    f.samples = 10000;
    f.real = 10000;
    f.value = 1;
    f.lower = 0.9996;
    f.upper = 1;
    const Json::Value v = dummy_fraction_json(f);
    EXPECT_EQ((std::vector<std::string>{ "interval", "samples", "source", "value" }),
              v.getMemberNames());
    EXPECT_EQ(1.0, v["value"].asDouble());
    ASSERT_EQ(2u, v["interval"].size());
    EXPECT_EQ(0.9996, v["interval"][0].asDouble());
    EXPECT_EQ(1.0, v["interval"][1].asDouble());
    EXPECT_EQ(10000u, v["samples"].asUInt64());
    EXPECT_EQ("sampled", v["source"].asString());
}

// --pattern-max-checked-entries (unchecked candidates tested one by one), as `metagraph
// pattern` and the server read it (the same Config): default 50, any integer in [0, 1000],
// refused at start-up beyond it or when not an integer; through the real loader the
// capabilities state it (the counts it makes exact: PatternRoute.UnmaskedTinyBlocksAreExact)
TEST(PatternMaskUnmasked, CheckedEntriesFlag) {
    const std::string dir = make_dir("checked_flag");
    const std::string built = build_masked(dir, "basic", "stat");
    const std::string stripped = strip(dir, built);
    const std::string anno = dir + "/anno.column.annodbg";
    {
        auto graph = load(built);
        annot::ColumnCompressed<> annotation(graph->max_index());
        std::vector<uint64_t> rows(graph->max_index());
        std::iota(rows.begin(), rows.end(), 0);
        annotation.add_labels(rows, { "records" });
        annotation.serialize(anno);
    }
    std::ofstream(dir + "/request.json") << "{}";
    auto config_of = [&](const std::string &graph, std::vector<std::string> flags) {
        std::vector<std::string> args = { "pattern" };
        args.insert(args.end(), flags.begin(), flags.end());
        for (const std::string &a : { std::string("-i"), graph, std::string("-a"), anno,
                                      dir + "/request.json" }) {
            args.push_back(a);
        }
        return make_config(args);
    };
    EXPECT_EQ(50u, config_of(stripped, {})->pattern_max_checked_entries);
    EXPECT_EQ(50u, pattern_limits(*config_of(stripped, {})).max_checked_entries);
    for (const char *value : { "0", "7", "1000" }) {
        EXPECT_EQ(std::stoull(value), pattern_limits(*config_of(
                stripped, { "--pattern-max-checked-entries", value })).max_checked_entries);
    }
    for (const char *value : { "1001", "-1", "x", "5.5", "" }) {
        EXPECT_DEATH(config_of(stripped, { "--pattern-max-checked-entries", value }),
                     "--pattern-max-checked-entries must be an integer in \\[0, 1000\\]")
            << value;
    }

    // through the loader: the capabilities state the limit
    auto absent = initialize_annotated_dbg(*config_of(stripped, {}));
    const PatternLimits on = pattern_limits(*config_of(stripped, {}));
    const PatternLimits off = pattern_limits(*config_of(stripped,
                                                        { "--pattern-max-checked-entries", "0" }));
    EXPECT_EQ(50u, pattern_capabilities_json(absent.get(), on, false)["caps"]
                           ["max_checked_entries"].asUInt64());
    EXPECT_EQ(0u, pattern_capabilities_json(absent.get(), off, false)["caps"]
                          ["max_checked_entries"].asUInt64());
}

} // namespace
