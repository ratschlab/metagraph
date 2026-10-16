#include "gtest/gtest.h"
#include <algorithm>
#include <map>
#include <random>
#include <set>
#include "graph/traversal/resolve.hpp"
#include "tests/annotation/coordinate_test_graph.hpp"
#include "cli/resolve.hpp"

#if defined(_DNA_GRAPH)
namespace {
using namespace mtg::graph::traversal;
using mtg::test::coordinate_test_graph;
const size_t kK = 11;
std::vector<KmerInterval> iv(std::initializer_list<std::pair<uint64_t, uint64_t>> l) {
    std::vector<KmerInterval> out;
    for (auto [a,b] : l) out.push_back({a,b});
    return out;
}
// These three selection/identity tests are extracted from Gunnar's prototype.
TEST(Resolve, MaxSupportScalesWithManyLabels) {
    const size_t num_labels = 1000;
    const uint64_t num_kmers = 5000;
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = num_kmers;
    profile.graph_runs = iv({ {0, num_kmers} });
    std::mt19937 gen(7);
    for (size_t i = 0; i < num_labels; ++i) {
        LabelProfile lp;
        lp.label.column = static_cast<uint64_t>(i);
        lp.label.name = "L" + std::to_string(i);
        // every label covers [0, 1000); each ends somewhere in [1000, 5000)
        uint64_t end = 1000 + gen() % (num_kmers - 1000);
        lp.runs = iv({ {0, end} });
        lp.kmers_supported = end;
        profile.labels.push_back(lp);
    }
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        profile.candidates.push_back({ profile.labels[l].runs[0], { l } });
    }

    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.max_seeds = 1;
    policy.min_block_bp = kK;   // min_kmers == 1
    auto selection = select_seeds(profile, std::string(num_kmers + kK - 1, 'A'), policy, false);
    ASSERT_EQ(1u, selection.seeds.size());
    // all 1000 labels cover [0, b) for b = the smallest run end
    uint64_t smallest_end = num_kmers;
    for (const auto &lp : profile.labels) {
        smallest_end = std::min(smallest_end, lp.runs[0].end);
    }
    EXPECT_EQ((KmerInterval{ 0, smallest_end }), selection.seeds[0].kmers);
    EXPECT_EQ(num_labels, selection.seeds[0].labels.size());

    // a length floor trades labels for length: the interval must be at least 3000
    // k-mers, so only the labels reaching that far can support it
    policy.min_block_bp = 3000 + kK - 1;
    auto longer = select_seeds(profile, std::string(num_kmers + kK - 1, 'A'), policy, false);
    ASSERT_EQ(1u, longer.seeds.size());
    EXPECT_GE(longer.seeds[0].kmers.size(), 3000u);
    EXPECT_LT(longer.seeds[0].labels.size(), num_labels);
    for (const auto &name : longer.seeds[0].labels) {
        const auto &lp = *std::find_if(profile.labels.begin(), profile.labels.end(),
                                       [&](const LabelProfile &x) { return x.label.name == name; });
        EXPECT_GE(lp.runs[0].end, longer.seeds[0].kmers.end);
    }
}

// Later max_support rounds must see the runs CLIPPED by what earlier rounds took.
// With A on [0,10) and B on [0,5), the first round takes the breadth-optimal [0,5)
// (2 labels); the second must then return [5,10), not a stub at the original run's end.
TEST(Resolve, MaxSupportClipsRunsBetweenRounds) {
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = 10;
    profile.graph_runs = iv({ {0, 10} });
    for (auto [name, end] : { std::pair<const char*, uint64_t>{ "A", 10 },
                              std::pair<const char*, uint64_t>{ "B", 5 } }) {
        LabelProfile lp;
        lp.label.column = static_cast<uint64_t>(profile.labels.size());
        lp.label.name = name;
        lp.runs = iv({ {0, end} });
        lp.kmers_supported = end;
        profile.labels.push_back(lp);
    }

    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.max_seeds = 2;
    policy.min_block_bp = kK;          // min_kmers == 1
    policy.merge_overlapping = false;  // keep the rounds' intervals as chosen
    auto sel = select_seeds(profile, std::string(10 + kK - 1, 'A'), policy, false);
    ASSERT_EQ(2u, sel.seeds.size());
    EXPECT_EQ((KmerInterval{ 0, 5 }), sel.seeds[0].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A", "B" }), sel.seeds[0].labels);
    EXPECT_EQ((KmerInterval{ 5, 10 }), sel.seeds[1].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A" }), sel.seeds[1].labels);
    // the two seeds tile the supported region without overlapping
    EXPECT_TRUE(sel.seeds[0].overlaps_with.empty());
    EXPECT_TRUE(sel.seeds[1].overlaps_with.empty());
}

// seed ids: canonical orientation and label order do not change the id; release does
TEST(Resolve, SeedId) {
    std::string s = "ACGTTGCAAGT";
    std::string rc = "ACTTGCAACGT";
    EXPECT_EQ(make_seed_id("r1", s, true, { "b", "a" }), make_seed_id("r1", rc, true, { "a", "b" }));
    EXPECT_NE(make_seed_id("r1", s, false, { "a" }), make_seed_id("r1", rc, false, { "a" }));
    EXPECT_NE(make_seed_id("r1", s, true, { "a" }), make_seed_id("r2", s, true, { "a" }));
    // length-prefixed labels: {"a,b"} != {"a","b"}
    EXPECT_NE(make_seed_id("", s, false, { "a,b" }), make_seed_id("", s, false, { "a", "b" }));
    EXPECT_EQ("x3o2x5", encode_runs(iv({ {0, 3}, {5, 10} }), 10));
}



TEST(Resolve, LiteralMembershipAcrossBackendsAndRecordBoundaries) {
    std::mt19937 rng(670);
    for (int trial = 0; trial < 12; ++trial) {
        std::string query(75, 'A');
        for (char &c : query) c = "ACGT"[rng() % 4];
        std::vector<mtg::test::CoordinateRecord> records{
            {query.substr(0, 32), "A", "same"}, {query.substr(43), "A", "same"},
            {query.substr(18, 41), "B", "same"}, {"AAAAAAAAAAAA", "C", "other"}};
        for (bool diff : {false,true}) for (bool brwt : {false,true}) {
            auto graph = coordinate_test_graph(5, records, diff, brwt);
            ResolveOptions options; options.discover = true;
            auto profile = resolve_support(*graph, query, options);
            std::map<std::string,std::vector<bool>> expected;
            std::vector<bool> present(query.size()-4);
            for (auto &record : records) {
                auto &mask=expected[record.column]; mask.resize(present.size());
                for (size_t i=0;i<mask.size();++i) {
                    bool match=record.sequence.find(query.substr(i,5))!=std::string::npos;
                    mask[i]=mask[i]||match; present[i]=present[i]||match;
                }
            }
            EXPECT_EQ(profile.graph_runs, runs_of(present.size(),[&](auto i){return present[i];}));
            size_t nonempty=0;
            for (const auto &[name,mask] : expected) nonempty+=std::count(mask.begin(),mask.end(),true)>0;
            ASSERT_EQ(profile.labels.size(),nonempty);
            for (const auto &label : profile.labels) {
                const auto &mask=expected.at(label.label.name);
                EXPECT_EQ(label.kmers_supported,std::count(mask.begin(),mask.end(),true));
                EXPECT_EQ(label.runs,runs_of(mask.size(),[&](auto i){return mask[i];}));
            }
            SelectionPolicy policy;policy.max_seeds=1;policy.min_block_bp=9;
            policy.merge_overlapping=false;
            auto chosen=select_seeds(profile,query,policy,false);
            size_t best_count=0;KmerInterval best;
            for (size_t a=0;a<present.size();++a) for(size_t b=a+5;b<=present.size();++b) {
                size_t count=0;
                for(const auto &[name,mask]:expected)
                    count+=std::all_of(mask.begin()+a,mask.begin()+b,[](bool x){return x;});
                if(count>best_count || (count==best_count && count && b-a>best.size())) {
                    best_count=count;best={a,b};
                }
            }
            if(!best_count) EXPECT_TRUE(chosen.seeds.empty());
            else {ASSERT_EQ(chosen.seeds.size(),1u);EXPECT_EQ(chosen.seeds[0].kmers,best);}
        }
    }
    // Membership may span different records. It must not be called a record trace.
    auto graph=coordinate_test_graph(3,{{"AACG","F","one"},{"ACGTT","F","two"}});
    ResolveOptions options;options.labels={"F"};
    EXPECT_EQ(resolve_support(*graph,"AACGTT",options).labels[0].runs,iv({{0,4}}));
}

TEST(Resolve, DiscoveryLimitsAndExplicitLabels) {
    auto graph=coordinate_test_graph(3,{{"ACGTTGCA","A","a"},{"ACGTTGCA","B","b"},
                                      {"ACGTTGCA","C","c"},{"ACG","D","d"}});
    ResolveOptions o;o.discover=true;o.discover_max_labels=2;
    auto p=resolve_support(*graph,"ACGTTGCA",o);
    ASSERT_TRUE(p.labels_truncated);EXPECT_EQ(p.labels_truncated->total,4u);
    EXPECT_EQ(p.labels_truncated->dropped_full_length,1u);
    EXPECT_EQ(p.labels[0].label.name,"A");EXPECT_EQ(p.labels[1].label.name,"B");
    o.discover=false;o.labels={"D","B"};p=resolve_support(*graph,"ACGTTGCA",o);
    EXPECT_EQ(p.labels[0].kmers_supported,1u);EXPECT_EQ(p.labels[1].kmers_supported,6u);
    o.labels={"A","A"};EXPECT_THROW(resolve_support(*graph,"ACGTTGCA",o),std::invalid_argument);
}

TEST(Resolve, RepeatedSequenceAndGapsAcrossRowBatches) {
    std::mt19937 rng(67017);
    std::string sequence(1200, 'A');
    for (char &c : sequence) c = "ACGT"[rng() % 4];
    const std::string query = sequence + "NNNNNNNNN" + sequence + sequence.substr(500, 500);
    const std::vector<mtg::test::CoordinateRecord> records{
        {sequence.substr(0, 700), "A", "a"}, {sequence.substr(500), "B", "b"}};
    for (bool diff : {false, true}) for (bool brwt : {false, true}) {
        auto graph = coordinate_test_graph(7, records, diff, brwt);
        for (bool discover : {false, true}) {
            ResolveOptions options;
            options.discover = discover;
            if (!discover) options.labels = {"B", "A"};
            const auto profile = resolve_support(*graph, query, options);
            ASSERT_EQ(profile.labels.size(), records.size());
            for (const auto &label : profile.labels) {
                const auto &record = *std::find_if(records.begin(), records.end(), [&](const auto &r) {
                    return r.column == label.label.name;
                });
                std::vector<bool> expected(query.size() - 6);
                for (size_t i = 0; i < expected.size(); ++i)
                    expected[i] = record.sequence.find(query.substr(i, 7)) != std::string::npos;
                EXPECT_EQ(label.kmers_supported, std::count(expected.begin(), expected.end(), true));
                EXPECT_EQ(label.runs, runs_of(expected.size(), [&](size_t i) { return expected[i]; }));
            }
        }
    }
}

TEST(Resolve, RequestRejectsUnsupportedSemantics) {
    auto graph=coordinate_test_graph(3,{{"ACGTTGCA","A","a"}});
    Json::Value request;request["sequence"]="ACGTTGCA";request["discover"]=Json::objectValue;
    auto result=mtg::cli::process_resolve_request(request,*graph);
    EXPECT_EQ(result["support"].asString(),"kmer");
    EXPECT_EQ(result["labels"][0]["kind"].asString(),"column");
    request["support"]="trace";
    EXPECT_THROW(mtg::cli::process_resolve_request(request,*graph),std::invalid_argument);
    request.removeMember("support");request["discover"]["kind"]="header";
    EXPECT_THROW(mtg::cli::process_resolve_request(request,*graph),std::invalid_argument);
    request["discover"].removeMember("kind");request["bounds"]["time_budget_ms"]=1;
    EXPECT_THROW(mtg::cli::process_resolve_request(request,*graph),std::invalid_argument);
}
} // namespace
#endif
